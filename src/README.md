# ディレクトリ構成

## Makefile

## 共通のライブラリ

行列用のファイル
- matrix.c
- matrix.h

誤差関数 α ダイバージェンスのファイル
- errfunc.c
- errfunc.h

反復アルゴリズムの停止条件のファイル
- stopkkt.c
- stopkkt.h

## 情報型更新

情報型更新ソルバー
- mur.c
- mur.h

情報型更新の数値実験ファイル

- test_mur.c

## 中津らのアルゴリズム

中津らのアルゴリズムのソルバー1

- nakatsu1.c
- nakatsu1.h

中津らのアルゴリズムのソルバー2

- nakatsu2.c
- nakatsu2.h

中津らのアルゴリズムの数値実験ファイル

- test_nakatsu1.c
- test_nakatsu2.c

## 提案法

提案法のソルバー

- mono.c
- mono.h

提案法のハイパーパラメータの数値実験ファイル

- test_beta.c

提案法の数値実験ファイル

- test_mono.c

## 数値実験用データセット

人工データ

- synth_1.csv
- synth_2.csv
- synth_3.csv

肺がん患者診断データ

- wdbc.csv

顔画像データ

- orl.csv

文書データ

- cluto_tr23.csv
