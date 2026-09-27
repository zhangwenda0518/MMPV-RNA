#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Title: 顶刊级批量正选择映射与 AlphaFold 模型智能渲染流水线 (终极完全体 + 标签自适应散开引擎)
@Description:
1. 从 CSV 自动解析正选择位点矩阵。
2. 彻底解决 PyMOL 无头模式命令行 argparse 传参的死锁与参数污染问题。
3. 智能模糊匹配 AlphaFold 输出长名文件。
4. 应用顶级学术期刊渲染管线，输出全部 7 大生化图景风格。
5. 【超大白底缓冲】使用 zoom buffer=25.0 绝对物理留白，拒绝切边。
6. 【高级四象限悬浮标签】采用动态向量数组 (Offset Matrix)，让聚集的突变点文字自动向四个方向散开，消除重叠。
7. 实时提取三维模型中的对应氨基酸名称，生成配套的注释图表 (.csv)。
8. 自动生成纯享版 .pml 命令行脚本，支持 GUI 拖拽微调。
"""

import os
import sys
import csv
import glob
import argparse
import traceback

try:
    import pymol
    from pymol import cmd, stored
except ImportError:
    print("❌ 致命错误: 必须在激活相应的 PyMOL Conda 环境中执行此脚本！", flush=True)
    sys.exit(1)

# =====================================================================
# 0. 底层通道隔离模块 
# =====================================================================
def get_clean_args():
    user_args =[]
    if '--' in sys.argv:
        user_args = sys.argv[sys.argv.index('--') + 1:]
    else:
        user_args = sys.argv[1:]
    return user_args

# =====================================================================
# 1. 组学生化数据解析模块
# =====================================================================
def parse_positive_sites(csv_path):
    positive_sites = {}
    if not os.path.exists(csv_path):
        print(f"\n❌[致命错误] 找不到输入的进化矩阵文件: {csv_path}", flush=True)
        cmd.quit(1)
        
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene = row.get('gene', '').strip()
                site = row.get('site', '').strip()
                fel = row.get('fel_selection', '').strip()
                meme = row.get('meme_marker', '').strip()
                
                if fel == 'positive' or (meme != '-' and meme != ''):
                    if gene not in positive_sites:
                        positive_sites[gene] = []
                    positive_sites[gene].append(site)
                    
        return positive_sites
    except Exception as e:
        print(f"❌ 解析 CSV 文件时底层异常: {e}", flush=True)
        cmd.quit(1)

# =====================================================================
# 2. 空间拓扑模式识别与顶点光照引擎配置
# =====================================================================
def find_structure_file(gene_name, search_dir="./"):
    core_id = gene_name.split('_')[-1].lower()
    patterns =[
        f"*_prot_{core_id}_model_0.cif",  
        f"*_{core_id}_unrelaxed_*.pdb",   
        f"{gene_name}.cif",               
        f"{gene_name}.pdb",               
        f"*{core_id}*model*.cif"          
    ]
    real_search_dir = os.path.expanduser(search_dir)
    for pattern in patterns:
        search_path = os.path.join(real_search_dir, pattern)
        matches = glob.glob(search_path)
        if matches:
            return matches[0] 
    return None

def apply_top_journal_settings():
    """ 出版级三维结构全局光影材质与字体基础设定 """
    cmd.bg_color("white")                
    cmd.set("orthoscopic", "on")         
    cmd.set("antialias", 2)              
    cmd.set("ray_trace_mode", 1)         
    cmd.set("ray_shadows", 0)            
    cmd.set("depth_cue", 1)              
    cmd.set("specular", 0.1)             
    cmd.set("ambient", 0.3)              
    cmd.set("cartoon_sampling", 14)      
    cmd.set("sphere_scale", 1.2)         
    
    # ⭐ 字体消肿与引线核心激活
    cmd.set("label_font_id", 5)          # 更改为常规无衬线体(Arial)，消除粗体的视觉黏连粘稠感
    cmd.set("label_size", 14)            # 字体大幅缩减至 14，留出更多呼吸空间
    cmd.set("label_shadow_mode", 2)      # 强阴影依然开启，保证白背景下的清晰度
    cmd.set("label_connector", 1)        # 开启引线 (必须打开，因为我们要把字体拉远)
    cmd.set("label_connector_ext_length", 1.0) 
    cmd.set("label_connector_width", 2.0)

# =====================================================================
# 3. 核心：四象限自适应散开偏置矩阵算法
# =====================================================================
def apply_staggered_labels(safe_obj, sites):
    """
    遍历该基因下的每一个突变位点，依次赋予不同方向的三维偏置量。
    避免空间共定位引发的标签文字挤压塌缩。
    """
    # X/Y 轴不同方向推远，Z 轴统一拔高 5.0 产生出屏脱离感
    offsets =[
        [3.5, 3.5, 5.0],   # 第一象限：右上[-3.5, 3.5, 5.0],  # 第二象限：左上[3.5, -3.5, 5.0],  # 第四象限：右下[-3.5, -3.5, 5.0], # 第三象限：左下[5.0, 0.0, 5.0],   # 正右[-5.0, 0.0, 5.0]   # 正左
    ]
    
    for idx, site in enumerate(sites):
        # 轮询调用散开序列向量
        offset = offsets[idx % len(offsets)]
        # 精确到该具体残基打入偏置量参数
        target_obj = f"{safe_obj} and resi {site} and name CA"
        cmd.alter(target_obj, f"label_position={offset}")

# =====================================================================
# 4. 产生离线手动 PML 指令文件
# =====================================================================
def write_manual_pml_script(gene, safe_obj, file_to_load, sites, site_str, output_dir):
    pml_content = f"""# =========================================================================
# 🧬 PyMOL Manual Rendering Script for Gene / Object: {gene}
# =========================================================================
# 【专家提示：极端重叠情况的手动标签位移指南】
# 1. 切换 PyMOL 鼠标模式为 "3-Button Editing" (点击右下角 Mouse Mode)
# 2. 按住 Ctrl + Shift (Mac为Cmd+Shift) 并用鼠标左键拖拽标签文字，任意排布！
#--------------------------------------------------------------------------

reinitialize
bg_color white
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0
set depth_cue, 1
set specular, 0.1
set ambient, 0.3
set cartoon_sampling, 14
set sphere_scale, 1.2

# 高阶引线标签字库载入 (轻量化字形)
set label_font_id, 5
set label_size, 14
set label_shadow_mode, 2
set label_color, black
set label_connector, 1
set label_connector_ext_length, 1.0
set label_connector_width, 2.0

# 挂载物理模型与建立视网膜焦距（加入超级 buffer 防止标签被切去边缘）
load {os.path.abspath(file_to_load)}, {safe_obj}
zoom {safe_obj}, buffer=25.0

# 全局选择集合定义
select pos_sites, {safe_obj} and resi {site_str}

# 生成悬浮引线标签
label pos_sites and name CA, "%s-%s" % (resn, resi)

#[自适应偏置矩阵引擎植入] 动态拉散靠近的突变热点标签
"""
    # 动态写入各偏移量参数至 PML 脚本
    offsets = [[3.5, 3.5, 5.0],[-3.5, 3.5, 5.0],[3.5, -3.5, 5.0], [-3.5, -3.5, 5.0],[5.0, 0.0, 5.0], [-5.0, 0.0, 5.0]]
    for idx, site in enumerate(sites):
        offset = offsets[idx % len(offsets)]
        pml_content += f"alter {safe_obj} and resi {site} and name CA, label_position={offset}\n"

    pml_content += f"""
# ---[Style A：常规半透表面透视 (默认开启)] ---
hide everything
show cartoon, {safe_obj}
color gray70, {safe_obj}
show surface, {safe_obj}
set transparency, 0.6
show spheres, pos_sites
color red, pos_sites
show labels
color black, pos_sites

# 其他风格以 # 注释隐去，需要在界面使用时自行删减井号：
# ---[Style D：经典小麦 SCI 基础沉浸款] -------------
# hide everything
# show cartoon, {safe_obj}
# color wheat, {safe_obj}
# show spheres, pos_sites
# color firebrick, pos_sites
# show labels

# ---[Style E：AlphaFold pLDDT 柔性区域热点上下文网络] ---
# hide everything
# show cartoon, {safe_obj}
# color orange, ({safe_obj}) and b < 50
# color yellow, ({safe_obj}) and b > 49.9 and b < 70
# color cyan, ({safe_obj}) and b > 69.9 and b < 90
# color blue, ({safe_obj}) and b > 89.9
# show spheres, pos_sites
# color magenta, pos_sites
# show labels
"""
    pml_file_path = os.path.join(output_dir, f"{gene}_ManualReplay_Settings.pml")
    try:
        with open(pml_file_path, "w", encoding="utf-8") as f:
            f.write(pml_content)
    except Exception:
        pass

# =====================================================================
# 5. 超线程批处理渲染主引擎 (Main Pipeline)
# =====================================================================
def main():
    print("\n[系统初始化] 🧬 PyMOL Top-Tier 结构进化渲染引擎全域启动...", flush=True)

    clean_args = get_clean_args()
    parser = argparse.ArgumentParser(description="PyMOL 全自动化 3D 正选位点图谱构建器")
    parser.add_argument('--input', required=True, help='输入的结合位点 CSV 结果文件矩阵')
    parser.add_argument('--output', default='./PyMOL_Figures', help='高清图谱集输出标的保存目录')
    parser.add_argument('--pdb_dir', default='./', help='存放 AlphaFold Cif/Pdb 元数据的索引目录')
    
    args = parser.parse_args(clean_args)

    abs_output = os.path.abspath(os.path.expanduser(args.output))
    os.makedirs(abs_output, exist_ok=True)

    val_input = os.path.expanduser(args.input)
    positive_sites_data = parse_positive_sites(val_input)
    
    if not positive_sites_data:
        print("⚠️ 警告: 捕获为 0。不包含命中的正选择位点矩阵，渲染管线正常休眠。", flush=True)
        cmd.quit(0)

    total = len(positive_sites_data)
    print(f"📊 数据就绪 | 锁定并瞄准 {total} 个演化基因宏 | 渲染流导向: {abs_output}\n", flush=True)

    annotation_records =[]

    for i, (gene, sites) in enumerate(positive_sites_data.items(), 1):
        file_to_load = find_structure_file(gene, args.pdb_dir)
        
        if not file_to_load:
            print(f"[{i}/{total}] ⚠️ 资产库中无法匹配识别 {gene}。跳跃。", flush=True)
            continue
            
        print(f"🚀 [{i}/{total}] -->[核心调度] 挂载资产: {os.path.basename(file_to_load)}")
        site_str = "+".join(sites)
        safe_obj = gene.replace('.', '_').replace('-', '_')
        label_target = f"{safe_obj} and resi {site_str} and name CA"
        
        try:
            # -------------------------------------------------------------
            #[基石信息收集模块] 
            # -------------------------------------------------------------
            cmd.reinitialize()
            cmd.load(file_to_load, safe_obj)
            stored.aa_dict = {}
            cmd.iterate(label_target, "stored.aa_dict[str(resi)] = resn")
            
            for site in sites:
                aa_name = stored.aa_dict.get(str(site), "Unknown/Missing")
                annotation_records.append({
                    "Gene_Object": gene,
                    "Structure_File": os.path.basename(file_to_load),
                    "Selection_Site": site,
                    "Amino_Acid": aa_name,
                    "Color_Style": "Dynamic"
                })

            write_manual_pml_script(gene, safe_obj, file_to_load, sites, site_str, abs_output)

            # -------------------------------------------------------------
            #[Style A：常规半透表面透视] 
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("gray70", safe_obj)
            cmd.show("surface", safe_obj)
            cmd.set("transparency", 0.6)         
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.show("spheres", "pos_sites")
            cmd.color("red", "pos_sites")        
            
            # 🔥 注入浮动引线引擎打签
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)")
            # 激活四象限散开算法
            apply_staggered_labels(safe_obj, sites)
            
            out_A = os.path.join(abs_output, f"{gene}_StyleA_Translucent.png")
            cmd.ray(1200, 1200)                  
            cmd.png(out_A)
            print(f"  ✅[Style A 图传] -> {os.path.basename(out_A)}", flush=True)

            # -------------------------------------------------------------
            #[Style B：抗原表位区实体块]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("surface", safe_obj)
            cmd.color("palecyan", safe_obj)          
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.color("red", "pos_sites")     
            
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)  
            
            out_B = os.path.join(abs_output, f"{gene}_StyleB_Surface.png")
            cmd.ray(1200, 1200)
            cmd.png(out_B)
            print(f"  ✅[Style B 图传] -> {os.path.basename(out_B)}", flush=True)

            # -------------------------------------------------------------
            #[Style C：生化侧链二级结构]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("cyan", f"{safe_obj} and ss h")     
            cmd.color("magenta", f"{safe_obj} and ss s")  
            cmd.color("lightpink", f"{safe_obj} and ss l")
            
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.show("sticks", "pos_sites")
            cmd.util.cbaw("pos_sites")  
            cmd.color("red", f"pos_sites and elem C") 
            
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)
            
            out_C = os.path.join(abs_output, f"{gene}_StyleC_Secondary.png")
            cmd.ray(1200, 1200)
            cmd.png(out_C)
            print(f"  ✅[Style C 图传] -> {os.path.basename(out_C)}", flush=True)

            # -------------------------------------------------------------
            #[Style D：经典小麦 SCI 基础沉浸款]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("wheat", safe_obj)             
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.show("spheres", "pos_sites")
            cmd.color("firebrick", "pos_sites") 
            
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)
            
            out_D = os.path.join(abs_output, f"{gene}_StyleD_Classic.png")
            cmd.ray(1200, 1200)
            cmd.png(out_D)
            print(f"  ✅[Style D 图传] -> {os.path.basename(out_D)}", flush=True)
            
            # -------------------------------------------------------------
            # ⭐[Style E：AlphaFold pLDDT 柔性区域热点上下文网络]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("orange", f"({safe_obj}) and b < 50")               
            cmd.color("yellow", f"({safe_obj}) and b > 49.9 and b < 70")    
            cmd.color("cyan", f"({safe_obj}) and b > 69.9 and b < 90")     
            cmd.color("blue", f"({safe_obj}) and b > 89.9")                 
            
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.show("spheres", "pos_sites")
            cmd.color("magenta", "pos_sites") 
            
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)
            
            out_E = os.path.join(abs_output, f"{gene}_StyleE_AF_Context.png")
            cmd.ray(1200, 1200)
            cmd.png(out_E)
            print(f"  ✅[Style E 图传] -> {os.path.basename(out_E)}", flush=True)

            # -------------------------------------------------------------
            # ⭐[Style F：近邻 5Å 生态微环境物理相互网]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("gray80", safe_obj)
            
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.select("neighbors", f"byres ({safe_obj} within 5.0 of pos_sites) and not pos_sites")
            
            cmd.show("sticks", "neighbors")
            cmd.color("gray40", "neighbors") 
            cmd.util.cbaw("neighbors")      
            
            cmd.show("sticks", "pos_sites")
            cmd.util.cbaw("pos_sites")
            cmd.color("red", f"pos_sites and elem C") 
            
            cmd.set("label_color", "black")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)
            
            out_F = os.path.join(abs_output, f"{gene}_StyleF_Microenvironment.png")
            cmd.ray(1200, 1200)
            cmd.png(out_F)
            print(f"  ✅[Style F 图传] -> {os.path.basename(out_F)}", flush=True)

            # -------------------------------------------------------------
            # ⭐[Style G：宇宙暗场聚类荧光构象足迹图]
            # -------------------------------------------------------------
            cmd.reinitialize()
            apply_top_journal_settings()
            cmd.bg_color("black")                
            cmd.set("ray_shadows", 0)            
            cmd.load(file_to_load, safe_obj)
            cmd.zoom("all", buffer=25.0)
            
            cmd.hide("everything")
            cmd.show("cartoon", safe_obj)
            cmd.color("gray20", safe_obj)            
            
            cmd.select("pos_sites", f"{safe_obj} and resi {site_str}")
            cmd.show("mesh", "pos_sites")        
            cmd.set("mesh_width", 2.0)           
            cmd.color("green", "pos_sites")      
            
            cmd.set("label_color", "yellow")
            cmd.label(label_target, "'%s-%s' % (resn, resi)") 
            apply_staggered_labels(safe_obj, sites)
            
            out_G = os.path.join(abs_output, f"{gene}_StyleG_MeshFootprint.png")
            cmd.ray(1200, 1200)
            cmd.png(out_G)
            print(f"  ✅[Style G 图传] -> {os.path.basename(out_G)}", flush=True)

        except Exception as e:
            print(f"  ❌ 截断警告！渲染节点脱靶: {gene} -> 异常流: {e}", flush=True)
            traceback.print_exc()

    # --- 阶段 D：将注释列表外置为数据字典 CSV ---
    annotation_csv_path = os.path.join(abs_output, "site_annotation_mapping.csv")
    if annotation_records:
        try:
            with open(annotation_csv_path, 'w', newline='', encoding='utf-8-sig') as csvfile:
                fieldnames = annotation_records[0].keys()
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for row in annotation_records:
                    writer.writerow(row)
            print(f"\n📝 附属数据模块激活！已生成靶点全览表")
        except Exception:
             pass

    cmd.quit(0)

main()
