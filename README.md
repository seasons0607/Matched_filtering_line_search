# Matched_filtering_line_search
- このスクリプトは、Matched filtering line search法 ([Miyazaki et al. 2016, PASJ, 68, 100](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true)) を用いて、輝線/吸収線の検出の有意性を検証するスクリプトです。
- 本スクリプトを用いて論文を執筆する場合は、[Miyazaki et al. 2016, PASJ, 68, 100](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true)とInoue et al., 2024, MNRAS, ??, ?? を引用してください。
- バグ等を見つけたら、お気軽にご連絡ください。
--  📧 inoue *at* cr.scphys.kyoto-u.ac.jp

## How to use
1. 必要なファイル
   1. 観測されたスペクトルデータをテキストファイル化したもの。
   2. 観測されたスペクトルデータをフィットする際に用いるrmf、arfファイル。
   3. 観測されたスペクトルデータのcontinuumのみをフィットし、テキストファイル化したもの。
   4. 観測されたスペクトルデータのcontinuumのみをフィットした際の、xcmファイル。


   を用意します。
2. main 関数
3. スクリプトを実行します。
   
>[!WARNING]
>本スクリプトの`calc_fwhm`関数は、NICERのrmfファイルを前提としています。別の衛星のrmfファイルに用いる際には、その衛星のrmfファイルの構造用に書き換えが必要な場合があります。また、FWHMのエネルギー依存性も、NICER/すざくの場合には3次の多項式で近似できますが、別の衛星ではより高次の式が必要な可能性があります。

## Example 
使用例として、NICERにより観測されたRS CVn型連星UX Ariの観測データに、本手法を用います。
1. まず必要なファイル一式を揃えます。
   1. 観測されたスペクトルをテキストファイル化
      iplot
      wd observed_spectrum.txt
3.
4.
5.


## Enviroment
Python 3 と 6.32 以降のHeasoftを前提としています。

## History
2024-10-04 new version 0.0 is created.

### Table of Content
[How to use](#how-to-use)
