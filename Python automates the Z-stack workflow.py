#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------------------------------------
功能（在原有脚本基础上）：
1) 读取 Z-stack（单文件或文件夹序列），分离 R/G/B；
2) 生成每通道的 XY/XZ/YZ 最大投影；可选分位数自动对比度；
3) 生成 RGB 合成图（支持 --red-scale 降低红色显示强度，红色=死亡细胞不突出）；
4) 针对绿色通道的 XY 投影，进行取向分析（结构张量），导出：
   - 0–180° 半圆玫瑰图（论文风格）
   - angles.csv（角度分布）
   - anisotropy.csv（各向异性指数 0–1）
5) 批处理：对目录中的多个样本（文件或子文件夹）逐个运行，并输出汇总表。

依赖：
  numpy, tifffile, imageio, matplotlib, scikit-image, pandas
"""

import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff


# 基础 I/O 
def is_image_file(p: Path) -> bool:
    return p.suffix.lower() in {".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp"}


def read_any(path: Path):
    suf = path.suffix.lower()
    if suf in (".tif", ".tiff"):
        from tifffile import imread
        return imread(str(path))
    else:
        import imageio.v2 as imageio
        return imageio.imread(str(path))


def stack_from_folder(folder: Path):
    files = sorted([f for f in folder.iterdir() if f.is_file() and is_image_file(f)])
    if not files:
        raise FileNotFoundError(f"No image files found in: {folder}")
    arrs = [read_any(f) for f in files]
    return np.stack(arrs, axis=0)  # (Z, Y, X, 3) or (Z, Y, X)


def load_input(input_path: str):
    """标准化输入形状：
       - ('rgbz', (Z,Y,X,3))  | ('grayz', (Z,Y,X))
       - ('rgb2d', (Y,X,3))   | ('cz', (C,Z,Y,X)) | ('zc', (Z,C,Y,X))
    """
    p = Path(input_path)
    if p.is_dir():
        vol = stack_from_folder(p)
        if vol.ndim == 4:
            return 'rgbz', vol
        return 'grayz', vol
    else:
        arr = read_any(p)
        if arr.ndim == 2:
            return 'grayz', arr[None, ...]
        if arr.ndim == 3:
            if arr.shape[-1] in (3, 4):
                return 'rgb2d', arr[..., :3]
            return 'grayz', arr
        if arr.ndim == 4:
            if arr.shape[-1] in (3, 4):
                return 'rgbz', arr[..., :3]
            if arr.shape[0] <= 6 and arr.shape[1] > 6:
                return 'cz', arr
            if arr.shape[1] <= 6 and arr.shape[0] > 6:
                return 'zc', arr
        raise ValueError(f"Unsupported input shape: {arr.shape}")


#  投影 and 对比度 
def mip(a, axis):
    return np.max(a, axis=axis)


def orthos_from_rgb_vol(Rvol, Gvol, Bvol):
    R_XY, G_XY, B_XY = mip(Rvol, 0), mip(Gvol, 0), mip(Bvol, 0)   # (Y,X)
    R_XZ, G_XZ, B_XZ = mip(Rvol, 1), mip(Gvol, 1), mip(Bvol, 1)   # (Z,X)
    R_YZ, G_YZ, B_YZ = mip(Rvol, 2), mip(Gvol, 2), mip(Bvol, 2)   # (Z,Y)
    return (R_XY, G_XY, B_XY), (R_XZ, G_XZ, B_XZ), (R_YZ, G_YZ, B_YZ)


def to01(a):
    a = a.astype(np.float32, copy=False)
    lo, hi = float(np.min(a)), float(np.max(a))
    if hi == lo:
        return np.zeros_like(a, dtype=np.float32)
    return (a - lo) / (hi - lo)


def auto_contrast(img, p_low=0.5, p_high=99.5):
    img = img.astype(np.float32, copy=False)
    lo, hi = np.percentile(img, [p_low, p_high])
    if hi <= lo:
        lo, hi = float(np.min(img)), float(np.max(img))
        if hi == lo:
            return np.zeros_like(img, dtype=np.float32)
    out = (img - lo) / (hi - lo + 1e-12)
    return np.clip(out, 0.0, 1.0)


def compose_rgb(R, G, B, scale_r=0.5, scale_g=1.0, scale_b=1.0):
    R = np.clip(R * scale_r, 0, 1)
    G = np.clip(G * scale_g, 0, 1)
    B = np.clip(B * scale_b, 0, 1)
    return (np.stack([R, G, B], axis=-1) * 255).astype(np.uint8)


#  Save
def save_all(XY_rgb8, XZ_rgb8, YZ_rgb8, chans, out_dir, prefix="mip"):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    import imageio.v2 as imageio
    for k, v in chans.items():
        t = (out / f"{prefix}_{k}.tif").as_posix()
        p = (out / f"{prefix}_{k}.png").as_posix()
        tiff.imwrite(t, v.astype(np.float32))
        imageio.imwrite(p, (np.clip(v, 0, 1) * 255).astype(np.uint8))
    tiff.imwrite((out / f"{prefix}_XY_RGB.tif").as_posix(), XY_rgb8, photometric='rgb')
    imageio.imwrite((out / f"{prefix}_XY_RGB.png").as_posix(), XY_rgb8)
    tiff.imwrite((out / f"{prefix}_XZ_RGB.tif").as_posix(), XZ_rgb8, photometric='rgb')
    imageio.imwrite((out / f"{prefix}_XZ_RGB.png").as_posix(), XZ_rgb8)
    tiff.imwrite((out / f"{prefix}_YZ_RGB.tif").as_posix(), YZ_rgb8, photometric='rgb')
    imageio.imwrite((out / f"{prefix}_YZ_RGB.png").as_posix(), YZ_rgb8)


#  XY 取向（绿色） 
def xy_orientation_and_exports(G_XY, out_dir, sample_name,
                               sigma=2.0, rose_bins=18,
                               mask_thresh='percentile', mask_pct=60.0):
    """
    - 对 GREEN XY 投影进行结构张量分析；
    - 以掩膜（percentile 或 Otsu）去掉背景；
    - 导出 0–180° 半圆玫瑰图、角度 CSV、各向异性 CSV。
    """
    # 仅导入 structure_tensor；特征值采用后备解析函数，兼容旧版 scikit-image
    from skimage.feature import structure_tensor
    from skimage import filters
    import matplotlib.pyplot as plt
    import pandas as pd

    def st_eigs(Axx, Axy, Ayy):
        # 2x2 结构张量的解析特征值（新版 API 缺失时使用）
        tmp = np.sqrt(((Axx - Ayy) * 0.5) ** 2 + Axy ** 2)
        l1 = 0.5 * (Axx + Ayy) + tmp
        l2 = 0.5 * (Axx + Ayy) - tmp
        return l1, l2

    outp = Path(out_dir) / "xy_orientation"
    outp.mkdir(parents=True, exist_ok=True)

    img = to01(G_XY)

    # 背景掩膜
    if mask_thresh == 'percentile':
        thr = np.percentile(img, mask_pct)
    else:
        try:
            thr = filters.threshold_otsu(img)
        except Exception:
            thr = np.percentile(img, 60.0)
    mask = img > thr

    # 结构张量 and 各向异性
    Axx, Axy, Ayy = structure_tensor(img, sigma=sigma)
    l1, l2 = st_eigs(Axx, Axy, Ayy)
    with np.errstate(divide='ignore', invalid='ignore'):
        anis_map = (l1 - l2) / (l1 + l2 + 1e-12)
    anis_map[~mask] = np.nan
    anis_value = np.nanmean(anis_map)

    # 角度（0–180°）
    theta = 0.5 * np.arctan2(2 * Axy, (Axx - Ayy))
    theta_deg = np.mod(np.degrees(theta), 180.0)
    theta_deg = theta_deg[mask & np.isfinite(theta_deg)]

    #  0–180° 半圆玫瑰图
    counts, bin_edges = np.histogram(theta_deg, bins=rose_bins, range=(0, 180))
    theta_c = np.deg2rad((bin_edges[:-1] + bin_edges[1:]) / 2)
    width = np.deg2rad(180.0 / rose_bins)
    rmax = (counts.max() * 1.10) if counts.size else 1

    fig = plt.figure(figsize=(4.2, 3.6), dpi=300)
    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location('E')     # 0° 在右侧
    ax.set_theta_direction(1)           # 逆时针增
    ax.set_thetamin(0); ax.set_thetamax(180)

    bars = ax.bar(theta_c, counts, width=width, bottom=0,
                  align='center', facecolor='black', edgecolor='black',
                  linewidth=0.6, alpha=1.0)

    ax.set_thetagrids(range(0, 181, 30))
    ax.set_ylim(0, rmax)
    ax.grid(alpha=0.25)
    ax.set_title(f"XY Orientation - {sample_name}", pad=10)
    ax.set_yticklabels([]) 

    fig.savefig((outp / f"{sample_name}_xy_rose.png").as_posix(),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    # CSV 导出
    import pandas as pd
    pd.DataFrame({"angles_deg": theta_deg}).to_csv(
        (outp / f"{sample_name}_xy_angles.csv").as_posix(), index=False)
    pd.DataFrame({"anisotropy": [anis_value]}).to_csv(
        (outp / f"{sample_name}_xy_anisotropy.csv").as_posix(), index=False)

    return anis_value


#  单样本流程 
def process_single_sample(input_path: Path, out_dir: Path,
                          auto_contrast_flag=False, p_low=0.5, p_high=99.5,
                          red_scale=0.5, do_xy_orient=False,
                          orient_sigma=2.0, rose_bins=18,
                          mask_thresh='percentile', mask_pct=60.0):
    mode, data = load_input(str(input_path))

    # 统一为 (Z,Y,X) 的 R/G/B
    if mode == 'rgbz':
        Rvol, Gvol, Bvol = data[..., 0], data[..., 1], data[..., 2]
    elif mode == 'rgb2d':
        Rvol = data[..., 0][None, ...]
        Gvol = data[..., 1][None, ...]
        Bvol = data[..., 2][None, ...]
    elif mode == 'cz':
        Rvol = data[0]
        Gvol = data[1] if data.shape[0] > 1 else np.zeros_like(data[0])
        Bvol = data[2] if data.shape[0] > 2 else np.zeros_like(data[0])
    elif mode == 'zc':
        Rvol = data[:, 0]
        Gvol = data[:, 1] if data.shape[1] > 1 else np.zeros_like(data[:, 0])
        Bvol = data[:, 2] if data.shape[1] > 2 else np.zeros_like(data[:, 0])
    elif mode == 'grayz':
        Rvol = Gvol = Bvol = data
    else:
        raise ValueError(f"Unhandled mode {mode}")

    # 三向投影
    (R_XY, G_XY, B_XY), (R_XZ, G_XZ, B_XZ), (R_YZ, G_YZ, B_YZ) = orthos_from_rgb_vol(
        Rvol.astype(np.float32), Gvol.astype(np.float32), Bvol.astype(np.float32)
    )

    # 归一化 + 可选自动对比度
    def proc(ch):
        chn = to01(ch)
        if auto_contrast_flag:
            chn = auto_contrast(chn, p_low, p_high)
        return chn

    R_XY, G_XY, B_XY = proc(R_XY), proc(G_XY), proc(B_XY)
    R_XZ, G_XZ, B_XZ = proc(R_XZ), proc(G_XZ), proc(B_XZ)
    R_YZ, G_YZ, B_YZ = proc(R_YZ), proc(G_YZ), proc(B_YZ)

    # RGB 合成（降低红色显示权重）
    XY_rgb8 = compose_rgb(R_XY, G_XY, B_XY, scale_r=red_scale)
    XZ_rgb8 = compose_rgb(R_XZ, G_XZ, B_XZ, scale_r=red_scale)
    YZ_rgb8 = compose_rgb(R_YZ, G_YZ, B_YZ, scale_r=red_scale)

    # 保存
    chans = {
        'R_XY': R_XY, 'G_XY': G_XY, 'B_XY': B_XY,
        'R_XZ': R_XZ, 'G_XZ': G_XZ, 'B_XZ': B_XZ,
        'R_YZ': R_YZ, 'G_YZ': G_YZ, 'B_YZ': B_YZ,
    }
    pref = "mip_autocontrast" if auto_contrast_flag else "mip_minmax"
    save_all(XY_rgb8, XZ_rgb8, YZ_rgb8, chans, out_dir, prefix=pref)

    # XY 取向分析（绿色）
    anis_value = None
    if do_xy_orient:
        sample_name = input_path.stem if input_path.is_file() else input_path.name
        anis_value = xy_orientation_and_exports(
            G_XY, out_dir, sample_name,
            sigma=orient_sigma, rose_bins=rose_bins,
            mask_thresh=mask_thresh, mask_pct=mask_pct
        )

    return anis_value


#  批处理 
def iter_samples_in_dir(top: Path):
    for p in sorted(top.iterdir()):
        if p.name.startswith('.') or p.name.startswith('_'):
            continue
        if p.is_file() and is_image_file(p):
            yield p.stem, p
        elif p.is_dir():
            imgs = [f for f in p.iterdir() if f.is_file() and is_image_file(f)]
            if imgs:
                yield p.name, p


def main():
    ap = argparse.ArgumentParser(description="Orthoview + XY orientation (GREEN) + red down-weighted composites")
    ap.add_argument("--input", required=True, help="单个栈文件/文件夹；或在 --batch 下为样本总目录")
    ap.add_argument("--save-dir", required=True, help="输出目录")
    # 对比度/显示
    ap.add_argument("--auto-contrast", action="store_true", help="启用分位数自动对比度")
    ap.add_argument("--p-low", type=float, default=0.5)
    ap.add_argument("--p-high", type=float, default=99.5)
    ap.add_argument("--red-scale", type=float, default=0.5, help="RGB 合成时红色缩放（仅显示）")
    # 取向
    ap.add_argument("--xy-orientation", action="store_true", help="对绿色 XY 投影做取向分析")
    ap.add_argument("--orient-sigma", type=float, default=2.0, help="结构张量高斯 σ（像素）")
    ap.add_argument("--rose-bins", type=int, default=18, help="半圆玫瑰图的分箱数")
    ap.add_argument("--mask-thresh", choices=["percentile","otsu"], default="percentile")
    ap.add_argument("--mask-pct", type=float, default=60.0, help="百分位掩膜阈值（保留高于该百分位的像素）")
    # 批处理
    ap.add_argument("--batch", action="store_true", help="把 --input 目录视为多个样本集合")
    ap.add_argument("--summary-csv", type=str, default="", help="与 --xy-orientation 一起使用，输出各向异性汇总")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_root = Path(args.save_dir); out_root.mkdir(parents=True, exist_ok=True)

    if not args.batch:
        anis = process_single_sample(
            in_path, out_root,
            auto_contrast_flag=args.auto_contrast,
            p_low=args.p_low, p_high=args.p_high,
            red_scale=args.red_scale,
            do_xy_orient=args.xy_orientation,
            orient_sigma=args.orient_sigma, rose_bins=args.rose_bins,
            mask_thresh=args.mask_thresh, mask_pct=args.mask_pct
        )
        if args.xy_orientation and args.summary_csv:
            import pandas as pd
            pd.DataFrame([{"sample": in_path.stem if in_path.is_file() else in_path.name,
                           "anisotropy": anis}]).to_csv(args.summary_csv, index=False)
    else:
        if not in_path.is_dir():
            raise ValueError("--batch 需要 --input 指向一个包含样本的目录")
        summary = []
        for name, spath in iter_samples_in_dir(in_path):
            sample_out = out_root / name
            sample_out.mkdir(parents=True, exist_ok=True)
            try:
                anis = process_single_sample(
                    spath, sample_out,
                    auto_contrast_flag=args.auto_contrast,
                    p_low=args.p_low, p_high=args.p_high,
                    red_scale=args.red_scale,
                    do_xy_orient=args.xy_orientation,
                    orient_sigma=args.orient_sigma, rose_bins=args.rose_bins,
                    mask_thresh=args.mask_thresh, mask_pct=args.mask_pct
                )
                if args.xy_orientation:
                    summary.append({"sample": name, "anisotropy": anis})
                print(f"[OK] {name}")
            except Exception as e:
                print(f"[FAIL] {name}: {e}")
        if args.xy_orientation and args.summary_csv:
            import pandas as pd
            pd.DataFrame(summary).to_csv(args.summary_csv, index=False)


if __name__ == "__main__":
    main()