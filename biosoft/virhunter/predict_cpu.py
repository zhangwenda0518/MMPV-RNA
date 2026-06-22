#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Grigorii Sukhorukov, Macha Nikolski
"""
VirHunter CPU prediction script.
Predicts viral contigs from a FASTA file using pretrained NN and RF models.
Runs on CPU only (no GPU or Ray required).
"""
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["TF_XLA_FLAGS"] = "--tf_xla_cpu_global_jit"
# loglevel : 0 all printed, 1 I not printed, 2 I and W not printed, 3 nothing printed
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import fire
import yaml
import tensorflow as tf
import numpy as np
from Bio import SeqIO
import pandas as pd
import ray
from utils import preprocess as pp
from pathlib import Path
from models import model_5, model_7, model_10
from joblib import load
import psutil
from typing import List, Union


def predict_nn(ds_path, nn_weights_path, length, n_cpus=1, batch_size=256):
    """
    Breaks down contigs into fragments
    and uses pretrained neural networks to give predictions for fragments
    """
    try:
        pid = psutil.Process(os.getpid())
        pid.cpu_affinity(range(n_cpus))
    except AttributeError:
        print("cpu allocation is not working properly. This will not impact the analysis results but may increase the runtime")

    try:
        seqs_ = list(SeqIO.parse(ds_path, "fasta"))
    except FileNotFoundError:
        raise Exception("test dataset was not found. Change ds variable")

    out_table = {
        "id": [],
        "length": [],
        "fragment": [],
        "pred_plant_5": [],
        "pred_vir_5": [],
        "pred_bact_5": [],
        "pred_plant_7": [],
        "pred_vir_7": [],
        "pred_bact_7": [],
        "pred_plant_10": [],
        "pred_vir_10": [],
        "pred_bact_10": [],
    }
    if not seqs_:
        raise ValueError("All sequences were smaller than length of the model")
    test_fragments = []
    test_fragments_rc = []
    ray.init(num_cpus=n_cpus, num_gpus=0, include_dashboard=False)
    for seq in seqs_:
        fragments_, fragments_rc, _ = pp.fragmenting([seq], length, max_gap=0.8,
                                                     sl_wind_step=int(length / 2))
        test_fragments.extend(fragments_)
        test_fragments_rc.extend(fragments_rc)
        for j in range(len(fragments_)):
            out_table["id"].append(seq.id)
            out_table["length"].append(len(seq.seq))
            out_table["fragment"].append(j)
    it = pp.chunks(test_fragments, int(len(test_fragments) / n_cpus + 1))
    test_encoded = np.concatenate(ray.get([pp.one_hot_encode.remote(s) for s in it]))
    it = pp.chunks(test_fragments_rc, int(len(test_fragments_rc) / n_cpus + 1))
    test_encoded_rc = np.concatenate(ray.get([pp.one_hot_encode.remote(s) for s in it]))
    ray.shutdown()

    for model, s in zip([model_5.model(length), model_7.model(length), model_10.model(length)], [5, 7, 10]):
        model.load_weights(Path(nn_weights_path, f"model_{s}_{length}.h5"))
        prediction = model.predict([test_encoded, test_encoded_rc], batch_size)
        out_table[f"pred_plant_{s}"].extend(list(prediction[..., 0]))
        out_table[f"pred_vir_{s}"].extend(list(prediction[..., 1]))
        out_table[f"pred_bact_{s}"].extend(list(prediction[..., 2]))
    return pd.DataFrame(out_table).round(3)


def predict_rf(df, rf_weights_path, length):
    """
    Using predictions by predict_nn and weights of a trained RF classifier gives a single prediction for a fragment
    """
    clf = load(Path(rf_weights_path, f"RF_{length}.joblib"))
    X = df[
        ["pred_plant_5", "pred_vir_5", "pred_plant_7", "pred_vir_7", "pred_plant_10", "pred_vir_10", ]]
    y_pred = clf.predict(X)
    mapping = {0: "plant", 1: "virus", 2: "bacteria"}
    df["RF_decision"] = np.vectorize(mapping.get)(y_pred)
    prob_classes = clf.predict_proba(X)
    df["RF_pred_plant"] = prob_classes[..., 0]
    df["RF_pred_vir"] = prob_classes[..., 1]
    df["RF_pred_bact"] = prob_classes[..., 2]
    return df


def predict_contigs(df):
    """
    Based on predictions of predict_rf for fragments gives a final prediction for the whole contig
    """
    df = (
        df.groupby(["id", "length", 'RF_decision'], sort=False)
        .size()
        .unstack(fill_value=0)
    )
    df = df.reset_index()
    df = df.reindex(['id', 'length', 'virus', 'plant', 'bacteria'], axis=1).fillna(value=0)
    conditions = [
        (df['virus'] > df['plant']) & (df['virus'] > df['bacteria']),
        (df['plant'] > df['virus']) & (df['plant'] > df['bacteria']),
        (df['bacteria'] >= df['plant']) & (df['bacteria'] >= df['virus']),
    ]
    choices = ['virus', 'plant', 'bacteria']
    df['decision'] = np.select(conditions, choices, default='bacteria')
    df = df.loc[:, ['id', 'length', 'virus', 'plant', 'bacteria', 'decision']]
    df = df.rename(columns={'virus': '# viral fragments', 'bacteria': '# bacterial fragments', 'plant': '# plant fragments'})
    df['# viral / # total'] = (df['# viral fragments'] / (df['# viral fragments'] + df['# bacterial fragments'] + df['# plant fragments'])).round(3)
    df = df.sort_values(by='# viral fragments', ascending=False)
    return df


def predict(
    input: Union[str, List[str]],
    weights: str,
    out_dir: str,
    cpu: int = 2,
    return_viral: bool = True,
    length: int = 750,
    config: str = None
):
    """
    Predicts viral contigs from the fasta file

    Arguments:
    input: Path(s) to the input file(s) with sequences for prediction (fasta format).
           Can be a single string or list of strings.
    weights: Path to the folder with weights of pretrained NN and RF weights.
             This folder should contain two subfolders 500 and 1000.
             Each of them contains corresponding weight.
    out_dir: Path to the folder, where to store output. You should create it.
    cpu: Number of CPUs to use (default: 2).
    return_viral: Return contigs annotated as viral by virhunter (fasta format) (default: True).
    length: Do predictions only for contigs > l. We suggest default l=750, as it was tested in the paper (default: 750).
    config: Optional config file path (for backward compatibility).
    """
    # 如果提供了config文件，优先使用config文件（向后兼容）
    if config is not None:
        with open(config, "r") as yamlfile:
            cf = yaml.load(yamlfile, Loader=yaml.FullLoader)

        input = cf["predict"]["test_ds"]
        weights = cf["predict"]["weights"]
        out_dir = cf["predict"]["out_path"]
        cpu = cf["predict"]["n_cpus"]
        return_viral = cf["predict"]["return_viral"]
        length = cf["predict"]["limit"]

    # 处理input参数
    if isinstance(input, str):
        input = [input]
    elif isinstance(input, list):
        pass
    else:
        raise ValueError('input should be a string or list of strings')

    # 验证输入参数
    assert len(input) > 0, 'input cannot be empty'
    for ts in input:
        assert Path(ts).exists(), f'{ts} does not exist'
    assert Path(weights).exists(), f'{weights} does not exist'
    assert isinstance(length, int), 'length should be an integer'
    assert isinstance(cpu, int), 'cpu should be an integer'

    # 创建输出目录
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    for ts in input:
        dfs_fr = []
        dfs_cont = []
        for l_ in (500, 1000):
            print(f'starting prediction for {Path(ts).name} for fragment length {l_}')
            df = predict_nn(
                ds_path=ts,
                nn_weights_path=weights,
                length=l_,
                n_cpus=cpu,
            )
            df = predict_rf(
                df=df,
                rf_weights_path=weights,
                length=l_,
            )
            dfs_fr.append(df.round(3))
            df = predict_contigs(df).round(3)
            dfs_cont.append(df)
            print('prediction finished')

        # 根据长度过滤并合并结果
        df_500 = dfs_fr[0][(dfs_fr[0]['length'] >= length) & (dfs_fr[0]['length'] < 1500)]
        df_1000 = dfs_fr[1][(dfs_fr[1]['length'] >= 1500)]
        df = pd.concat([df_1000, df_500], ignore_index=True)
        pred_fr = Path(out_dir, f"{Path(ts).stem}_predicted_fragments.csv")
        df.to_csv(pred_fr)

        df_500 = dfs_cont[0][(dfs_cont[0]['length'] >= length) & (dfs_cont[0]['length'] < 1500)]
        df_1000 = dfs_cont[1][(dfs_cont[1]['length'] >= 1500)]
        df = pd.concat([df_1000, df_500], ignore_index=True)
        pred_contigs = Path(out_dir, f"{Path(ts).stem}_predicted.csv")
        df.to_csv(pred_contigs)

        if return_viral:
            viral_ids = list(df[df["decision"] == "virus"]["id"])
            seqs_ = list(SeqIO.parse(ts, "fasta"))
            viral_seqs = [s_ for s_ in seqs_ if s_.id in viral_ids]
            SeqIO.write(viral_seqs, Path(out_dir, f"{Path(ts).stem}_viral.fasta"), 'fasta')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
        print(__doc__)
        print("""
Usage: python predict_cpu.py [OPTIONS]

Options:
  --input PATH       Input FASTA file(s), string or list of strings.
  --weights PATH     Folder containing pretrained NN and RF weights.
  --out_dir PATH     Output directory path.
  --cpu INT          Number of CPUs to use (default: 2).
  --return_viral     Return contigs annotated as viral in a separate FASTA (default: True).
  --length INT       Only predict for contigs longer than this threshold (default: 750).
  --config PATH      Optional YAML config file path (for backward compatibility).

Examples:
  python predict_cpu.py --input=test.fasta --weights=./weights/tomato --out_dir=./output
  python predict_cpu.py --input=test.fasta --weights=./weights/tomato --out_dir=./output --cpu=4 --length=500
  python predict_cpu.py --config=configs/predict_config.yaml  (backward compatible)
        """)
        sys.exit(0)
    fire.Fire(predict)
