import streamlit as st
import psi4
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile


au2kcal=psi4.constants.hartree2kcalmol

st.set_page_config(page_title="Psi4計算", layout="wide")
st.title("簡易版量子化学計算webアプリ")

# 入力選択
input_type = st.radio("入力形式を選んでください：", ("SMILES", "XYZファイル"))

if input_type == "SMILES":
    smiles = st.text_input("SMILESを入力してください：")
else:
    xyz_file = st.file_uploader("XYZファイルをアップロードしてください：", type=["xyz"])

# NMR核種の選択
calc_type = st.radio("実行したい計算を選択してください：", ("エネルギー計算", "構造最適化", "振動数計算", "分子軌道可視化"))

# 計算方法の選択
method = st.radio("計算方法を選択してください：", ("hf", "mp2", "b3lyp"))

#既定関数の選択
basis = st.radio("基底関数を選択してください：", ("3-21g","6-31g","6-31g(d)","sto-3g", "cc-pvdz", "cc-pvtz", "aug-cc-pvtz"))

# 計算ボタン
if st.button("計算実行！"):
    with st.spinner("計算中..."):

        # 構造の処理
        if input_type == "SMILES":
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            xyz = Chem.MolToXYZBlock(mol)
        else:
            xyz = xyz_file.read().decode("utf-8")

        # 一時ファイルとして保存
        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".xyz") as temp:
            temp.write(xyz)
            temp.flush()
            molname = temp.name

        # Psi4設定
        psi4.core.set_output_file("psi4_output.dat", False)
        psi4.set_memory('10GB')


        # Psi4計算
        mol_obj = psi4.geometry(xyz)
        method_for_calc = method + "/" + basis

        if calc_type=="エネルギー計算":
            energy=psi4.energy(method_for_calc, molecule=mol_obj)
            success_message="Energy: "+str(round(energy*au2kcal,2))+" kcal/mol"
            st.success(success_message)

        elif calc_type=="構造最適化":
            # 最適化を実行
            psi4.optimize(method_for_calc, molecule=mol_obj)
            st.success("最適化完了！")
            # 最適化後の分子構造をXYZ形式で取り出す
            # 構造データをxyzに保存
            # 最適化後のXYZ文字列を取得
            xyz_full = mol_obj.save_string_xyz()
            # 行で分割
            xyz_lines = xyz_full.strip().splitlines()
            # 上3行を除いた部分だけ抽出
            xyz_trimmed = "\n".join(xyz_lines[3:])
            st.download_button(label="ダウンロード", file_name="optimized_struture.xyz", data=xyz)

        elif calc_type=="振動数計算":
            if method=="b3lyp":
                st.write("B3LYPによる振動数計算では文字コードエラーが発生します")
            else:
                psi4.set_options({"NORMAL_MODES_WRITE":True})
                psi4.optimize(method_for_calc, molecule=mol_obj)
                _, wfn=psi4.frequency(method_for_calc, molecule=mol_obj,return_wfn=True)
                st.write(wfn.frequencies().to_array())
                st.write("計算終了")
                st.write("avogadroで読み込む場合、拡張子を「.molden」に変更してください")

        elif calc_type=="分子軌道可視化":
            psi4.set_output_file("fchk.log")
            psi4.optimize(method_for_calc, molecule=mol_obj)
            _,wfn=psi4.energy(method_for_calc, molecule=mol_obj,return_wfn=True)
            psi4.fchk(wfn,"MO.fchk")
            st.write("計算終了")
            st.write("avogadroなどで読み込んでください")