# :bulb:Matched_filtering_line_search
- このスクリプトは、Matched filtering line search法 ([Miyazaki et al. 2016, PASJ, 68, 100](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true)) を用いて、輝線/吸収線の検出の有意性を検証するスクリプトです。
- 本スクリプトを用いて論文を執筆する場合は、[Miyazaki et al. 2016, PASJ, 68, 100](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true)とInoue et al., 2024, MNRAS, ??, ?? を引用してください。
- バグ等を見つけたら、お気軽にご連絡ください。
  -  📧 inoue *at* cr.scphys.kyoto-u.ac.jp

## :white_check_mark: How to use 
使用例として、NICERにより観測されたRS CVn型連星σ Gemの観測データの6.4 keV付近に見られた吸収構造 (Figure 2 of Inoue et al. 2024) に対して、本手法を用いて統計的有意性を検証します。

![spec_S1_absorption_paper_main.pdf](https://github.com/user-attachments/files/17375949/spec_S1_absorption_paper_main.pdf)


1. まず、必要なファイル一式を揃えます。
   1. 観測されたスペクトルをテキストファイル化します。本調査では6.4 keV付近の吸収線の有無を調べたいので、その周辺のバンド (5-8 keV) のデータのみを用います。
      ```
      data observed_spectrum.pi
      response observed_spectrum.rmf
      arf observed_spectrum.arf
      back observed_spectrum_3C50_bkg.pi
      ignore **-5.0 8.0-**
      iplot
      wd observed_spectrum.txt
      ```
   
   2. 観測されたスペクトルのContinuumのみをフィットします。本例では、鉄輝線を除いた5.0-6.2、7.2-8.0 keVを`bremss`モデルでフィットしています。
      ```
      data observed_spectrum.pi
      response observed_spectrum.rmf
      arf observed_spectrum.arf
      back observed_spectrum_3C50_bkg.pi
      ignore **-5.0 8.0-**
      ignore 6.2-7.2
      model bremss
      6.0
      0.1
      renorm
      fit
      ```
      フィットを完了したら、xcmファイルとして本モデルを保存します。
      ```
      save all continuum.xcm
      ```
      その後、continuumのモデルをテキストファイル化します。この際、観測されたスペクトルとContinuumのモデルの間でエネルギーのビン数を揃えるために、フィットに用いなかった6.2-7.2 keVのモデルデータもテキストファイルには含めます。
      ```
      notice all
      ignore **-5.0 8.0-**
      iplot
      wd continuum.txt
      ```
      これで、下準備は完了です。この過程で生成したファイルは全て、観測データやrmfファイルと同じディレクトリに置いてください。

> [!NOTE]
> 本手法はContinuumのモデルフィットが十分な精度で可能であることが前提のため、エネルギーバンドは対象とするライン±数keVの狭帯域に絞った方が良い場合が多いです。
2. main関数の中のパラメータを埋めていきます。まずは、rmfファイルや、1.で作成したテキストファイルのpathを記入します。
      ```  
      nicer_rmf_path = "/Desktop/1100380108/observed_spectrum.rmf"
      nicer_arf_path = "/Desktop/UX_Ari/1100380108/observed_spectrum.arf"
      nicer_spectrum_txt_path = "/Desktop/1100380108/observed_spectrum.txt"
      nicer_continuum_txt_path = "/Desktop/1100380108/continuum.txt"
      nicer_continuum_xcm_path = "/Desktop/1100380108/continuum.xcm"
      ```
      Exposureは、観測されたスペクトルと同じ値に設定します。シミュレーションの試行回数は`N=10000`程度は最低でも必要です。
      ``` 
      exposure = 1113 #sec
      trial_number = 10000
      ```
      最後に、1.(i)で観測データをテキストファイル化する際に指定したエネルギーバンドの下限値と上限値を入力します。
      ``` 
      continuum_energy_low = "5.0" #keV
      continuum_energy_upp = "8.0" #keV
      ```
3. スクリプトを実行します。ここからは基本的に待つだけです。
      ```
      python matched_filtering_line_search.py
      ```
   すると、まずは`calc_fwhm`関数が実行されます。この際、進行状況がステータスバーとして表示されます。NICERの場合は、１分ほどで終了します。
      ```
      Calculating FWHM of 1100380108_block012.rmf:  17%|████████████▊                      | 75/450 [00:09<00:47,  7.87it/s]
      ```
   計算が完了すると、FWHMのエネルギー依存性が出力されます。

   ![1200040104_block0_FWHM.pdf](https://github.com/user-attachments/files/17377393/1200040104_block0_FWHM.pdf)

   その後、`plot_spectrum_with_MC`関数が実行されます。まずは、fakeスペクトルの作成が開始し、この際も進行状況がステータスバーとして表示されます。`N=10000`の場合は、数時間ほどかかります。
      ```
      Generating fake spectra:  12%|███████████▋                                       | 12/100 [00:40<04:59,  3.40s/it]
      ```
   最後に、生成された全てのfakeスペクトルにMatched filterをかけます。この工程は、`N=10000`の場合は、数十分ほどかかります。
      ```
      Matched-filtering fake spectra: 100%|█████████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:02<00:00, 42.72it/s
      ```
   全ての計算が完了すると、観測されたスペクトルにMatched filterをかけたものとcontinuumの±1σ、±3σ を示す領域を示す図 (c.f. Figure 2 of [Miyazaki et al. 2016](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true): Figure 3 of Inoue et al. 2024) が出力され、対象とする構造が統計的に有意かどうか判定できます (軸ラベルや表示範囲は手動で微調整しましょう)。この場合は、吸収構造の中心が3σ下限値を下回っているため、3σ以上の有意度での検出を主張できます。
   
   ![Matched_filtered_spectrum.pdf](https://github.com/user-attachments/files/17377195/Matched_filtered_spectrum.pdf)

> [!TIP]
> 軸ラベルや表示範囲は手動で微調整しましょう。

## :white_check_mark: Functions
- `gaussian`
  - ガウス関数を書き下した関数です。
- `resolution_sigma`
  - FWHMのエネルギー依存性をフィットする際に用いる関数です。NICER/すざくでは、3次の多項式を用いています。
- `matched＿filter`
  - [Miyazaki et al. 2016](https://academic.oup.com/pasj/article/68/6/100/2664382?login=true) の式(1)に従い、Matched filterをかける関数です。
- `extract_integers`
  - ある範囲内の整数のみを取り出す関数です。図を出力する際の軸ラベルの設定に用います。
- `calc_fwhm`
  - rmfファイルから、FWHMの依存性を計算する関数です。`resolution_sigma`を用いて3次の多項式でフィットを行い、多項式の係数を出力します。
- `plot_spectrum_with_MC`
  - fakeスペクトルを作成し、Matched filterをかけた上で、その統計的分布を調べます。最終的に、観測されたスペクトルにMatched filterをかけたものとcontinuumの±1σ、±3σ を示す領域を示す図を出力します。
   
>[!WARNING]
>本スクリプトの`calc_fwhm`関数は、NICERのrmfファイルを前提としています。別の衛星のrmfファイルに用いる際には、その衛星のrmfファイルの構造用に書き換えが必要な場合があります。また、FWHMのエネルギー依存性も、NICER/すざくの場合には3次の多項式で近似できますが、別の衛星ではより高次の式が必要な可能性があります。


## :white_check_mark: Enviroment
Python 3 と 6.32 以降のHeasoftを前提としています。

## :white_check_mark: History
2024-10-15 new version 0.0 is created.

