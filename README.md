
 
``` python
import pandas as pd
from glob import glob
import numpy as np
```
 

 
## HTS Analysis
 

 
``` python
filepaths = glob("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\*.csv")
```
 

 
``` python
LOD_percentages = []
above_LOD_count = []
combined_dict = {}
for filepath in filepaths:
    dataframe_csv = pd.read_csv(filepath, encoding = "ISO-8859-1", usecols=lambda x: x not in ["Time", "T° 340"])
    column_name = filepath[-8:-4]
    slopes = dataframe_csv.apply(lambda x: np.polyfit(dataframe_csv.index, x, 1)[0])
    reshaped_df = pd.DataFrame(slopes.values.reshape(16,24), index=list('ABCDEFGHIJKLMNOP'), columns=list(range(1,25)))
    neg_mean = (reshaped_df[1].mean() + reshaped_df[2].mean())/2
    neg_std = pd.concat([reshaped_df[1],reshaped_df[2]]).std()
    LOD = (-1)*abs(neg_mean + 3*neg_std)
    pos_mean = (reshaped_df[23].mean() + reshaped_df[24].mean())/2
    LOD_percentage = (LOD*100)/pos_mean
    percentage_activity = slopes.apply(lambda x : (((x - neg_mean)*100)/pos_mean))
    LOD_percentages.append(LOD_percentage)
    above_LOD_count.append(percentage_activity.where(percentage_activity > LOD_percentage, np.nan).count())
    combined_dict[column_name] = percentage_activity.to_list()
combined_df = pd.DataFrame(data=combined_dict, index=slopes.index)
combined_df_display = pd.DataFrame(combined_df)
combined_df_display.loc["LOD percentage"] = LOD_percentages
combined_df_display.loc["Above LOD Count"] = above_LOD_count
combined_df_display
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>1-R1</th>
      <th>1-R2</th>
      <th>2-R1</th>
      <th>2-R2</th>
      <th>3-R1</th>
      <th>3-R2</th>
      <th>4-R1</th>
      <th>4-R2</th>
      <th>5-R1</th>
      <th>5-R2</th>
      <th>6-R1</th>
      <th>6-R2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A1</th>
      <td>-0.454741</td>
      <td>0.672701</td>
      <td>0.528338</td>
      <td>0.916917</td>
      <td>0.090003</td>
      <td>-0.469922</td>
      <td>11.035921</td>
      <td>-0.863511</td>
      <td>0.394905</td>
      <td>1.924535</td>
      <td>-0.166397</td>
      <td>0.163383</td>
    </tr>
    <tr>
      <th>A2</th>
      <td>0.192002</td>
      <td>0.672701</td>
      <td>0.835735</td>
      <td>-0.120880</td>
      <td>4.006912</td>
      <td>0.610120</td>
      <td>-0.606369</td>
      <td>1.576953</td>
      <td>-0.683859</td>
      <td>0.046463</td>
      <td>0.381735</td>
      <td>-0.391128</td>
    </tr>
    <tr>
      <th>A3</th>
      <td>-1.020640</td>
      <td>0.771902</td>
      <td>0.118476</td>
      <td>0.633882</td>
      <td>-0.601217</td>
      <td>0.610120</td>
      <td>-0.683985</td>
      <td>0.474808</td>
      <td>-0.298587</td>
      <td>0.515981</td>
      <td>-0.009788</td>
      <td>-0.549559</td>
    </tr>
    <tr>
      <th>A4</th>
      <td>-0.616426</td>
      <td>0.573501</td>
      <td>0.425873</td>
      <td>-0.875641</td>
      <td>-1.868452</td>
      <td>-0.137602</td>
      <td>-0.994446</td>
      <td>-0.706062</td>
      <td>-0.144477</td>
      <td>0.046463</td>
      <td>-0.009788</td>
      <td>0.480246</td>
    </tr>
    <tr>
      <th>A5</th>
      <td>-0.131370</td>
      <td>0.573501</td>
      <td>1.040666</td>
      <td>0.162156</td>
      <td>3.315693</td>
      <td>2.687125</td>
      <td>-0.528754</td>
      <td>0.553533</td>
      <td>-0.529750</td>
      <td>-0.110043</td>
      <td>5.236627</td>
      <td>-0.866422</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>P22</th>
      <td>-0.858955</td>
      <td>0.375101</td>
      <td>0.835735</td>
      <td>0.445191</td>
      <td>-1.407639</td>
      <td>0.111639</td>
      <td>-0.916830</td>
      <td>-0.391163</td>
      <td>0.394905</td>
      <td>0.281222</td>
      <td>0.146821</td>
      <td>0.084167</td>
    </tr>
    <tr>
      <th>P23</th>
      <td>103.913296</td>
      <td>101.559303</td>
      <td>101.149536</td>
      <td>100.734123</td>
      <td>98.358354</td>
      <td>100.472519</td>
      <td>102.699556</td>
      <td>103.761563</td>
      <td>99.487105</td>
      <td>104.279461</td>
      <td>103.430725</td>
      <td>103.777602</td>
    </tr>
    <tr>
      <th>P24</th>
      <td>104.640881</td>
      <td>102.948106</td>
      <td>105.760487</td>
      <td>101.960611</td>
      <td>102.736077</td>
      <td>104.294208</td>
      <td>100.759174</td>
      <td>108.721216</td>
      <td>103.570999</td>
      <td>107.018316</td>
      <td>106.014780</td>
      <td>106.787801</td>
    </tr>
    <tr>
      <th>LOD percentage</th>
      <td>0.802278</td>
      <td>1.963182</td>
      <td>2.376378</td>
      <td>1.899767</td>
      <td>5.746865</td>
      <td>0.425204</td>
      <td>2.401031</td>
      <td>0.771107</td>
      <td>0.935854</td>
      <td>0.885130</td>
      <td>1.592567</td>
      <td>0.804021</td>
    </tr>
    <tr>
      <th>Above LOD Count</th>
      <td>70.000000</td>
      <td>64.000000</td>
      <td>56.000000</td>
      <td>56.000000</td>
      <td>40.000000</td>
      <td>98.000000</td>
      <td>48.000000</td>
      <td>66.000000</td>
      <td>59.000000</td>
      <td>73.000000</td>
      <td>58.000000</td>
      <td>68.000000</td>
    </tr>
  </tbody>
</table>
<p>386 rows × 12 columns</p>
</div>
```
 
 

 
``` python
LOD_avg_cutoff = np.mean(LOD_percentages)
LOD_avg_cutoff
```

 
    1.716948778506576
 
 

 
``` python
LOD_percentages
```

 
    [0.8022782448852165,
     1.9631820982061519,
     2.3763783935249236,
     1.8997665139734345,
     5.746864931541303,
     0.425204050601918,
     2.401031132613852,
     0.7711074764966598,
     0.9358540780383339,
     0.8851304543629532,
     1.5925673455582974,
     0.8040206222758678]
 
 

 
``` python
remove_list = [str(x) + '1' for x in list('ABCDEFGHIJKLMNOP')]+[str(x) + '2' for x in list('ABCDEFGHIJKLMNOP')]+[str(x) + '23' for x in list('ABCDEFGHIJKLMNOP')]+[str(x) + '24' for x in list('ABCDEFGHIJKLMNOP')]
frames = []
for i in range(1,7):
    i = str(i)
    tmp_df = pd.DataFrame(combined_df[[i+'-R1',i+'-R2']]).rename(columns={i+'-R1':'R1', i+'-R2':'R2'})
    tmp_df = tmp_df.drop(remove_list)
    index_list = ['Plate '+i+'-' + x for x in tmp_df.index]
    tmp_df.index = index_list
    tmp_df['Is_Considered'] = (tmp_df['R1'] > 1.59) & (tmp_df['R2'] > 0.8)
    frames.append(tmp_df)
final_df = pd.concat(frames)
final_df['Rate Average'] = final_df[['R1','R2']].mean(axis=1)
final_df[final_df['Is_Considered']]
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>R1</th>
      <th>R2</th>
      <th>Is_Considered</th>
      <th>Rate Average</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Plate 1-D4</th>
      <td>9.084708</td>
      <td>10.394321</td>
      <td>True</td>
      <td>9.739514</td>
    </tr>
    <tr>
      <th>Plate 1-D10</th>
      <td>4.072455</td>
      <td>4.045508</td>
      <td>True</td>
      <td>4.058982</td>
    </tr>
    <tr>
      <th>Plate 1-D12</th>
      <td>46.676604</td>
      <td>150.861802</td>
      <td>True</td>
      <td>98.769203</td>
    </tr>
    <tr>
      <th>Plate 1-E4</th>
      <td>2.051386</td>
      <td>1.069502</td>
      <td>True</td>
      <td>1.560444</td>
    </tr>
    <tr>
      <th>Plate 1-E16</th>
      <td>2.293914</td>
      <td>3.648707</td>
      <td>True</td>
      <td>2.971311</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Plate 6-N18</th>
      <td>1.791220</td>
      <td>1.510051</td>
      <td>True</td>
      <td>1.650635</td>
    </tr>
    <tr>
      <th>Plate 6-O13</th>
      <td>4.531885</td>
      <td>6.817507</td>
      <td>True</td>
      <td>5.674696</td>
    </tr>
    <tr>
      <th>Plate 6-P4</th>
      <td>3.357314</td>
      <td>3.015150</td>
      <td>True</td>
      <td>3.186232</td>
    </tr>
    <tr>
      <th>Plate 6-P15</th>
      <td>1.634611</td>
      <td>1.113972</td>
      <td>True</td>
      <td>1.374291</td>
    </tr>
    <tr>
      <th>Plate 6-P16</th>
      <td>2.730877</td>
      <td>2.143777</td>
      <td>True</td>
      <td>2.437327</td>
    </tr>
  </tbody>
</table>
<p>110 rows × 4 columns</p>
</div>
```
 
 

 
``` python
final_df.sort_values('Rate Average', ascending=False)
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>R1</th>
      <th>R2</th>
      <th>Is_Considered</th>
      <th>Rate Average</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Plate 1-D12</th>
      <td>46.676604</td>
      <td>150.861802</td>
      <td>True</td>
      <td>98.769203</td>
    </tr>
    <tr>
      <th>Plate 6-B3</th>
      <td>135.300739</td>
      <td>57.277948</td>
      <td>True</td>
      <td>96.289344</td>
    </tr>
    <tr>
      <th>Plate 2-L17</th>
      <td>39.772654</td>
      <td>123.660004</td>
      <td>True</td>
      <td>81.716329</td>
    </tr>
    <tr>
      <th>Plate 2-B10</th>
      <td>83.115594</td>
      <td>76.298720</td>
      <td>True</td>
      <td>79.707157</td>
    </tr>
    <tr>
      <th>Plate 3-C7</th>
      <td>64.258199</td>
      <td>66.575798</td>
      <td>True</td>
      <td>65.416998</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Plate 4-J10</th>
      <td>-1.537753</td>
      <td>-3.855048</td>
      <td>False</td>
      <td>-2.696400</td>
    </tr>
    <tr>
      <th>Plate 4-D20</th>
      <td>-4.176672</td>
      <td>-1.257134</td>
      <td>False</td>
      <td>-2.716903</td>
    </tr>
    <tr>
      <th>Plate 6-G5</th>
      <td>-3.533500</td>
      <td>-3.322111</td>
      <td>False</td>
      <td>-3.427805</td>
    </tr>
    <tr>
      <th>Plate 2-H10</th>
      <td>-7.054115</td>
      <td>-1.158677</td>
      <td>False</td>
      <td>-4.106396</td>
    </tr>
    <tr>
      <th>Plate 2-E8</th>
      <td>-1.213577</td>
      <td>-7.762840</td>
      <td>False</td>
      <td>-4.488208</td>
    </tr>
  </tbody>
</table>
<p>1920 rows × 4 columns</p>
</div>
```
 
 

 
``` python
compound_df = pd.read_excel('C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\All 1A Data_EPA ECHO plate location.xlsx', usecols=['Dest plate', 'Dest well', 'EPA_SAMPLE_ID','PREFERRED_NAME'])
def combine_columns(row):
    if row['Dest well'][1] == '0':
        dest_well = row['Dest well'][0]+row['Dest well'][2]
    else:
        dest_well = row['Dest well']
    return "Plate " + str(row['Dest plate']) + "-" + dest_well

compound_df['ID'] = compound_df.apply(combine_columns, axis=1)
compound_df = compound_df.drop(columns=['Dest plate', 'Dest well'])
compound_df = compound_df.set_index("ID")
compound_df
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>EPA_SAMPLE_ID</th>
      <th>PREFERRED_NAME</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Plate 1-A3</th>
      <td>EPAPLT0439A09</td>
      <td>Imazamox</td>
    </tr>
    <tr>
      <th>Plate 1-B3</th>
      <td>EPAPLT0439A10</td>
      <td>Dimethyl succinate</td>
    </tr>
    <tr>
      <th>Plate 1-C3</th>
      <td>EPAPLT0439B09</td>
      <td>Toluene 2,4-diisocyanate</td>
    </tr>
    <tr>
      <th>Plate 1-D3</th>
      <td>EPAPLT0439B10</td>
      <td>1-Methylnaphthalene</td>
    </tr>
    <tr>
      <th>Plate 1-E3</th>
      <td>EPAPLT0439C09</td>
      <td>MEHP</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Plate 6-L22</th>
      <td>EPAPLT0536G01</td>
      <td>1-(Methylamino)-4-(2-hydroxyethylamino)anthraq...</td>
    </tr>
    <tr>
      <th>Plate 6-M22</th>
      <td>EPAPLT0536G04</td>
      <td>Sodium thiomethoxide</td>
    </tr>
    <tr>
      <th>Plate 6-N22</th>
      <td>EPAPLT0536G05</td>
      <td>Tris(tribromoneopentyl)phosphate</td>
    </tr>
    <tr>
      <th>Plate 6-O22</th>
      <td>EPAPLT0513A02</td>
      <td>Nickel chloride</td>
    </tr>
    <tr>
      <th>Plate 6-P22</th>
      <td>EPAPLT0485H05</td>
      <td>Fluopicolide</td>
    </tr>
  </tbody>
</table>
<p>1920 rows × 2 columns</p>
</div>
```
 
 

 
``` python
final_df = compound_df.join(final_df).sort_values('Rate Average', ascending=False)
```
 

 
``` python
final_df.plot.scatter(x='R1',y='R2')
```

 
    <Axes: xlabel='R1', ylabel='R2'>
 

 
![](images/c829017acf76badc83c7ed7adb783b8fcc405142.png)
 
 

 
``` python
final_df
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>EPA_SAMPLE_ID</th>
      <th>PREFERRED_NAME</th>
      <th>R1</th>
      <th>R2</th>
      <th>Is_Considered</th>
      <th>Rate Average</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Plate 1-D12</th>
      <td>EPAPLT0443B12</td>
      <td>Dichlone</td>
      <td>46.676604</td>
      <td>150.861802</td>
      <td>True</td>
      <td>98.769203</td>
    </tr>
    <tr>
      <th>Plate 6-B3</th>
      <td>EPAPLT0481G11</td>
      <td>9-Phenanthrol</td>
      <td>135.300739</td>
      <td>57.277948</td>
      <td>True</td>
      <td>96.289344</td>
    </tr>
    <tr>
      <th>Plate 2-L17</th>
      <td>EPAPLT0454H03</td>
      <td>Methylene blue</td>
      <td>39.772654</td>
      <td>123.660004</td>
      <td>True</td>
      <td>81.716329</td>
    </tr>
    <tr>
      <th>Plate 2-B10</th>
      <td>EPAPLT0449G03</td>
      <td>Diquat dibromide monohydrate</td>
      <td>83.115594</td>
      <td>76.298720</td>
      <td>True</td>
      <td>79.707157</td>
    </tr>
    <tr>
      <th>Plate 3-C7</th>
      <td>EPAPLT0458H12</td>
      <td>2,5-Di-tert-butylbenzene-1,4-diol</td>
      <td>64.258199</td>
      <td>66.575798</td>
      <td>True</td>
      <td>65.416998</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Plate 4-J10</th>
      <td>EPAPLT0466C03</td>
      <td>Dioctyl sebacate</td>
      <td>-1.537753</td>
      <td>-3.855048</td>
      <td>False</td>
      <td>-2.696400</td>
    </tr>
    <tr>
      <th>Plate 4-D20</th>
      <td>EPAPLT0470B08</td>
      <td>Disodium 4,4'-bis(2-sulfostyryl)biphenyl</td>
      <td>-4.176672</td>
      <td>-1.257134</td>
      <td>False</td>
      <td>-2.716903</td>
    </tr>
    <tr>
      <th>Plate 6-G5</th>
      <td>EPAPLT0483C03</td>
      <td>FR150011</td>
      <td>-3.533500</td>
      <td>-3.322111</td>
      <td>False</td>
      <td>-3.427805</td>
    </tr>
    <tr>
      <th>Plate 2-H10</th>
      <td>EPAPLT0450A03</td>
      <td>Ametryn</td>
      <td>-7.054115</td>
      <td>-1.158677</td>
      <td>False</td>
      <td>-4.106396</td>
    </tr>
    <tr>
      <th>Plate 2-E8</th>
      <td>EPAPLT0448H12</td>
      <td>Tetrac</td>
      <td>-1.213577</td>
      <td>-7.762840</td>
      <td>False</td>
      <td>-4.488208</td>
    </tr>
  </tbody>
</table>
<p>1920 rows × 6 columns</p>
</div>
```
 
 

 
``` python
final_LOD_cutoff_df = final_df[final_df["Is_Considered"]]
final_LOD_cutoff_df
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>EPA_SAMPLE_ID</th>
      <th>PREFERRED_NAME</th>
      <th>R1</th>
      <th>R2</th>
      <th>Is_Considered</th>
      <th>Rate Average</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Plate 1-D12</th>
      <td>EPAPLT0443B12</td>
      <td>Dichlone</td>
      <td>46.676604</td>
      <td>150.861802</td>
      <td>True</td>
      <td>98.769203</td>
    </tr>
    <tr>
      <th>Plate 6-B3</th>
      <td>EPAPLT0481G11</td>
      <td>9-Phenanthrol</td>
      <td>135.300739</td>
      <td>57.277948</td>
      <td>True</td>
      <td>96.289344</td>
    </tr>
    <tr>
      <th>Plate 2-L17</th>
      <td>EPAPLT0454H03</td>
      <td>Methylene blue</td>
      <td>39.772654</td>
      <td>123.660004</td>
      <td>True</td>
      <td>81.716329</td>
    </tr>
    <tr>
      <th>Plate 2-B10</th>
      <td>EPAPLT0449G03</td>
      <td>Diquat dibromide monohydrate</td>
      <td>83.115594</td>
      <td>76.298720</td>
      <td>True</td>
      <td>79.707157</td>
    </tr>
    <tr>
      <th>Plate 3-C7</th>
      <td>EPAPLT0458H12</td>
      <td>2,5-Di-tert-butylbenzene-1,4-diol</td>
      <td>64.258199</td>
      <td>66.575798</td>
      <td>True</td>
      <td>65.416998</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>Plate 5-B13</th>
      <td>EPAPLT0475E10</td>
      <td>1-Dodecyl-2-pyrrolidinone</td>
      <td>1.704833</td>
      <td>1.298511</td>
      <td>True</td>
      <td>1.501672</td>
    </tr>
    <tr>
      <th>Plate 2-L3</th>
      <td>EPAPLT0447E08</td>
      <td>Tri-allate</td>
      <td>1.757925</td>
      <td>1.199953</td>
      <td>True</td>
      <td>1.478939</td>
    </tr>
    <tr>
      <th>Plate 3-A19</th>
      <td>EPAPLT0463B08</td>
      <td>Bromoxynil</td>
      <td>1.818051</td>
      <td>0.942441</td>
      <td>True</td>
      <td>1.380246</td>
    </tr>
    <tr>
      <th>Plate 6-P15</th>
      <td>EPAPLT0487H03</td>
      <td>4-Butyloxyaniline</td>
      <td>1.634611</td>
      <td>1.113972</td>
      <td>True</td>
      <td>1.374291</td>
    </tr>
    <tr>
      <th>Plate 2-P4</th>
      <td>EPAPLT0448A02</td>
      <td>Benodanil</td>
      <td>1.757925</td>
      <td>0.822572</td>
      <td>True</td>
      <td>1.290249</td>
    </tr>
  </tbody>
</table>
<p>110 rows × 6 columns</p>
</div>
```
 
 

 
``` python
hits_dict = {}
Rate_avg = final_LOD_cutoff_df["Rate Average"]
for i in list(range(100,9,-10))+[5,3,2,1,0]:
    #print(r"Hits >= " + str(i) + r"% of Positive control : ",Rate_avg.where(Rate_avg >= i).count())
    hits_dict[str(i)+r"% and above"] = [Rate_avg.where(Rate_avg >= i).count(), Rate_avg.where(Rate_avg >= i).count()*100/final_df["Rate Average"].count()]
hits_df = pd.DataFrame.from_dict(hits_dict, orient="index", columns=["Number of Hits", "Toxcast Phase II library %"])
hits_df.index.name = "Activity compared to PR"
hits_df
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Number of Hits</th>
      <th>Toxcast Phase II library %</th>
    </tr>
    <tr>
      <th>Activity compared to PR</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>100% and above</th>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>90% and above</th>
      <td>2</td>
      <td>0.104167</td>
    </tr>
    <tr>
      <th>80% and above</th>
      <td>3</td>
      <td>0.156250</td>
    </tr>
    <tr>
      <th>70% and above</th>
      <td>4</td>
      <td>0.208333</td>
    </tr>
    <tr>
      <th>60% and above</th>
      <td>6</td>
      <td>0.312500</td>
    </tr>
    <tr>
      <th>50% and above</th>
      <td>7</td>
      <td>0.364583</td>
    </tr>
    <tr>
      <th>40% and above</th>
      <td>7</td>
      <td>0.364583</td>
    </tr>
    <tr>
      <th>30% and above</th>
      <td>9</td>
      <td>0.468750</td>
    </tr>
    <tr>
      <th>20% and above</th>
      <td>9</td>
      <td>0.468750</td>
    </tr>
    <tr>
      <th>10% and above</th>
      <td>17</td>
      <td>0.885417</td>
    </tr>
    <tr>
      <th>5% and above</th>
      <td>44</td>
      <td>2.291667</td>
    </tr>
    <tr>
      <th>3% and above</th>
      <td>71</td>
      <td>3.697917</td>
    </tr>
    <tr>
      <th>2% and above</th>
      <td>89</td>
      <td>4.635417</td>
    </tr>
    <tr>
      <th>1% and above</th>
      <td>110</td>
      <td>5.729167</td>
    </tr>
    <tr>
      <th>0% and above</th>
      <td>110</td>
      <td>5.729167</td>
    </tr>
  </tbody>
</table>
</div>
```
 
 

 
``` python
final_LOD_cutoff_df.to_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_HTS_analysis_above_LOD.xlsx")
```
 

 
``` python
final_df.count()
```

 
    EPA_SAMPLE_ID     1920
    PREFERRED_NAME    1920
    R1                1920
    R2                1920
    Is_Considered     1920
    Rate Average      1920
    dtype: int64
 
 

 
``` python
final_LOD_cutoff_df.count()
```

 
    EPA_SAMPLE_ID     110
    PREFERRED_NAME    110
    R1                110
    R2                110
    Is_Considered     110
    Rate Average      110
    dtype: int64
 
 

 
``` python
toxcast_library = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\toxcast_library.xlsx")
toxcast_library = toxcast_library[["PREFERRED NAME","SMILES"]]
toxcast_library.rename(columns={"PREFERRED NAME":"PREFERRED_NAME"}, inplace=True)
```
 

 
``` python
toxcast_library
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>PREFERRED_NAME</th>
      <th>SMILES</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Fluthiacet-methyl</td>
      <td>COC(=O)CSC1=C(Cl)C=C(F)C(=C1)N=C1SC(=O)N2CCCCN12</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2-(4-Chloro-2-methylphenoxy)acetic acid</td>
      <td>CC1=C(OCC(O)=O)C=CC(Cl)=C1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Thiophanate-methyl</td>
      <td>COC(=O)NC(=S)NC1=C(NC(=S)NC(=O)OC)C=CC=C1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Fludioxonil</td>
      <td>FC1(F)OC2=CC=CC(=C2O1)C1=CNC=C1C#N</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Acephate</td>
      <td>COP(=O)(NC(C)=O)SC</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6549</th>
      <td>Diethyl ethoxymethylenemalonate</td>
      <td>CCOC=C(C(=O)OCC)C(=O)OCC</td>
    </tr>
    <tr>
      <th>6550</th>
      <td>Dapsone</td>
      <td>NC1=CC=C(C=C1)S(=O)(=O)C1=CC=C(N)C=C1</td>
    </tr>
    <tr>
      <th>6551</th>
      <td>1-(4-Chlorophenyl)ethanone</td>
      <td>CC(=O)C1=CC=C(Cl)C=C1</td>
    </tr>
    <tr>
      <th>6552</th>
      <td>Flufenpyr-ethyl</td>
      <td>CCOC(=O)COC1=CC(N2N=CC(=C(C)C2=O)C(F)(F)F)=C(F...</td>
    </tr>
    <tr>
      <th>6553</th>
      <td>Bromadiolone</td>
      <td>OC(CC(C1=CC=CC=C1)C1=C(O)C2=C(OC1=O)C=CC=C2)C1...</td>
    </tr>
  </tbody>
</table>
<p>6554 rows × 2 columns</p>
</div>
```
 
 

 
``` python
final_df_with_mols = final_df.merge(toxcast_library, how='left', on='PREFERRED_NAME').drop_duplicates()
final_df_with_mols = final_df_with_mols.set_index("EPA_SAMPLE_ID")
```
 

 
``` python
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import DataStructs
import os
from rdkit import RDConfig
PandasTools.RenderImagesInAllDataFrames(images=True)
```
 

 
``` python
final_df_with_mols["SMILES"].fillna("C", inplace=True)
PandasTools.AddMoleculeColumnToFrame(final_df_with_mols, smilesCol='SMILES', molCol='MOL_OBJ', includeFingerprints=True)
```
 

 
``` python
final_df_with_mols_above_lod = final_df_with_mols[final_df_with_mols["Is_Considered"]]
final_df_with_mols_above_lod[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]]
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>PREFERRED_NAME</th>
      <th>R1</th>
      <th>R2</th>
      <th>Rate Average</th>
      <th>MOL_OBJ</th>
    </tr>
    <tr>
      <th>EPA_SAMPLE_ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>EPAPLT0443B12</th>
      <td>Dichlone</td>
      <td>46.676604</td>
      <td>150.861802</td>
      <td>98.769203</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0481G11</th>
      <td>9-Phenanthrol</td>
      <td>135.300739</td>
      <td>57.277948</td>
      <td>96.289344</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0454H03</th>
      <td>Methylene blue</td>
      <td>39.772654</td>
      <td>123.660004</td>
      <td>81.716329</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0449G03</th>
      <td>Diquat dibromide monohydrate</td>
      <td>83.115594</td>
      <td>76.298720</td>
      <td>79.707157</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0458H12</th>
      <td>2,5-Di-tert-butylbenzene-1,4-diol</td>
      <td>64.258199</td>
      <td>66.575798</td>
      <td>65.416998</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>EPAPLT0475E10</th>
      <td>1-Dodecyl-2-pyrrolidinone</td>
      <td>1.704833</td>
      <td>1.298511</td>
      <td>1.501672</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0447E08</th>
      <td>Tri-allate</td>
      <td>1.757925</td>
      <td>1.199953</td>
      <td>1.478939</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0463B08</th>
      <td>Bromoxynil</td>
      <td>1.818051</td>
      <td>0.942441</td>
      <td>1.380246</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0487H03</th>
      <td>4-Butyloxyaniline</td>
      <td>1.634611</td>
      <td>1.113972</td>
      <td>1.374291</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
    <tr>
      <th>EPAPLT0448A02</th>
      <td>Benodanil</td>
      <td>1.757925</td>
      <td>0.822572</td>
      <td>1.290249</td>
      <td style="text-align: center;"><div style="width: 200px; height: 200px" data-content="rdkit/molecule"><img src="data:image/png;base64,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" alt="Mol"/></div></td>
    </tr>
  </tbody>
</table>
<p>110 rows × 5 columns</p>
</div>
```
 
 

 
``` python
final_df_with_mols_above_lod[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]].to_html("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\above_LOD_2Y3_screen_data_with_structures.html")
```
 

 
``` python
final_df_with_mols[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]].to_html("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_screen_data_with_structures.html")
```
 

 
``` python
PandasTools.SaveXlsxFromFrame(final_df_with_mols_above_lod, "C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\above_LOD_2Y3_screen_data_with_structures.xlsx", molCol='MOL_OBJ')
```
 

 
``` python
```
 

 
## HTS Followup
 

 
``` python
ABS_followup_data = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ABS_followup.xlsx")
ROS_followup_data = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ROS_followup.xlsx")
follow_up_picklist = pd.read_csv("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14 Derek_ABS_picklist.csv")
```
 

 
``` python
follow_up_slopes = ABS_followup_data.apply(lambda x: np.polyfit(ABS_followup_data.index, x, 1)[0])
follow_up_reshaped_df = pd.DataFrame(follow_up_slopes.values.reshape(16,24), index=list('ABCDEFGHIJKLMNOP'), columns=list(range(1,25)))

neg_mean = follow_up_reshaped_df[23].mean()
neg_std = follow_up_reshaped_df[23].std()
pos_mean = follow_up_reshaped_df[22].mean()
percentage_activity = follow_up_slopes.apply(lambda x : (((x - neg_mean)*100)/pos_mean))
follow_up_reshaped_df_percentage_activity = follow_up_reshaped_df.apply(lambda x : (((x - neg_mean)*100)/pos_mean))
follow_up_comp_list = follow_up_picklist[0:48]
follow_up_comp_list['ID'] = follow_up_comp_list['Source Plate Barcode'].apply(lambda x : "Plate " + x[12])+"-"+follow_up_comp_list['Source Well']
follow_up_comp_list = follow_up_comp_list.merge(compound_df, how='left', on='ID').drop_duplicates()[['EPA_SAMPLE_ID','PREFERRED_NAME']]
follow_up_dict_abs = {"No Protein R1":[1,2,3], "No Protein R2":[4,5,6], "Reductase only R1":[7,8,9], "Reductase only R2":[10,11,12], "2Y3-Reductase R1":[13,14,15], "2Y3-Reductase R2":[16,17,18], "2Y3-Reductase R3":[19,20,21]}
for label in follow_up_dict_abs.keys():
    follow_up_comp_list[label] = follow_up_reshaped_df_percentage_activity[follow_up_dict_abs[label]].to_numpy().flatten('F')
```

 
    C:\Users\thisi\AppData\Local\Temp\ipykernel_27828\1258482499.py:10: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead

    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      follow_up_comp_list['ID'] = follow_up_comp_list['Source Plate Barcode'].apply(lambda x : "Plate " + x[12])+"-"+follow_up_comp_list['Source Well']
 
 

 
``` python
follow_up_slopes_ros = ROS_followup_data.apply(lambda x: np.polyfit(ROS_followup_data.index, x, 1)[0])
follow_up_reshaped_df_ros = pd.DataFrame(follow_up_slopes_ros.values.reshape(16,24), index=list('ABCDEFGHIJKLMNOP'), columns=list(range(1,25)))

neg_mean_ros = follow_up_reshaped_df_ros[20].mean()
neg_std_ros = follow_up_reshaped_df_ros[20].std()
pos_mean_ros = follow_up_reshaped_df_ros[19].mean()
follow_up_reshaped_df_ros_percentage_activity = follow_up_reshaped_df_ros.apply(lambda x : (((x - neg_mean_ros)*100)/pos_mean_ros))
follow_up_dict_ros = {"Reductase only ROS R1":[1,2,3], "Reductase only ROS R2":[4,5,6], "Reductase only ROS R3":[7,8,9], "2Y3-Reductase ROS R1":[10,11,12], "2Y3-Reductase ROS R2":[13,14,15], "2Y3-Reductase ROS R3":[16,17,18]}
for label in follow_up_dict_ros.keys():
    follow_up_comp_list[label] = follow_up_reshaped_df_ros_percentage_activity[follow_up_dict_ros[label]].to_numpy().flatten('F')
```
 

 
``` python
follow_up_comp_list['Avg No Protein'] = follow_up_comp_list[['No Protein R1','No Protein R2']].mean(axis=1)
follow_up_comp_list['Avg Reductase only'] = follow_up_comp_list[['Reductase only R1','Reductase only R2']].mean(axis=1)
follow_up_comp_list['Avg 2Y3-Reductase'] = follow_up_comp_list[['2Y3-Reductase R1','2Y3-Reductase R2','2Y3-Reductase R3']].mean(axis=1)
follow_up_comp_list['Avg Reductase only ROS'] = follow_up_comp_list[['Reductase only ROS R1','Reductase only ROS R2','Reductase only ROS R3']].mean(axis=1)
follow_up_comp_list['Avg 2Y3-Reductase ROS'] = follow_up_comp_list[['2Y3-Reductase ROS R1','2Y3-Reductase ROS R2','2Y3-Reductase ROS R3']].mean(axis=1)
follow_up_comp_list.sort_values("Avg 2Y3-Reductase", ascending=False)
```

 
```{=html}
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>EPA_SAMPLE_ID</th>
      <th>PREFERRED_NAME</th>
      <th>No Protein R1</th>
      <th>No Protein R2</th>
      <th>Reductase only R1</th>
      <th>Reductase only R2</th>
      <th>2Y3-Reductase R1</th>
      <th>2Y3-Reductase R2</th>
      <th>2Y3-Reductase R3</th>
      <th>Reductase only ROS R1</th>
      <th>Reductase only ROS R2</th>
      <th>Reductase only ROS R3</th>
      <th>2Y3-Reductase ROS R1</th>
      <th>2Y3-Reductase ROS R2</th>
      <th>2Y3-Reductase ROS R3</th>
      <th>Avg No Protein</th>
      <th>Avg Reductase only</th>
      <th>Avg 2Y3-Reductase</th>
      <th>Avg Reductase only ROS</th>
      <th>Avg 2Y3-Reductase ROS</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>9</th>
      <td>EPAPLT0454H03</td>
      <td>Methylene blue</td>
      <td>8.964091</td>
      <td>12.019307</td>
      <td>97.310773</td>
      <td>94.255556</td>
      <td>166.307749</td>
      <td>153.917148</td>
      <td>150.437596</td>
      <td>29.817660</td>
      <td>33.750273</td>
      <td>44.972652</td>
      <td>50.271349</td>
      <td>10.284370</td>
      <td>43.217234</td>
      <td>10.491699</td>
      <td>95.783164</td>
      <td>156.887498</td>
      <td>36.180195</td>
      <td>34.590984</td>
    </tr>
    <tr>
      <th>41</th>
      <td>EPAPLT0481G11</td>
      <td>9-Phenanthrol</td>
      <td>-0.371294</td>
      <td>0.053042</td>
      <td>97.819976</td>
      <td>99.856787</td>
      <td>155.529624</td>
      <td>157.311834</td>
      <td>152.729009</td>
      <td>6.764996</td>
      <td>8.455525</td>
      <td>12.610127</td>
      <td>34.752637</td>
      <td>35.842089</td>
      <td>53.707050</td>
      <td>-0.159126</td>
      <td>98.838381</td>
      <td>155.190155</td>
      <td>9.276883</td>
      <td>41.433925</td>
    </tr>
    <tr>
      <th>10</th>
      <td>EPAPLT0449G03</td>
      <td>Diquat dibromide monohydrate</td>
      <td>-2.917308</td>
      <td>-3.341643</td>
      <td>41.807670</td>
      <td>42.401740</td>
      <td>84.495836</td>
      <td>86.108312</td>
      <td>86.787249</td>
      <td>9.514240</td>
      <td>7.413886</td>
      <td>3.076570</td>
      <td>38.970421</td>
      <td>37.500173</td>
      <td>42.103876</td>
      <td>-3.129475</td>
      <td>42.104705</td>
      <td>85.797132</td>
      <td>6.668232</td>
      <td>39.524823</td>
    </tr>
    <tr>
      <th>36</th>
      <td>EPAPLT0481D10</td>
      <td>1,4-Dihydroxy-2-naphthoic acid</td>
      <td>-4.020580</td>
      <td>-3.765979</td>
      <td>37.394579</td>
      <td>41.468201</td>
      <td>78.470270</td>
      <td>75.330186</td>
      <td>74.651249</td>
      <td>0.334157</td>
      <td>6.235639</td>
      <td>1.608030</td>
      <td>29.179016</td>
      <td>26.738849</td>
      <td>29.228536</td>
      <td>-3.893280</td>
      <td>39.431390</td>
      <td>76.150568</td>
      <td>2.725942</td>
      <td>28.382133</td>
    </tr>
    <tr>
      <th>11</th>
      <td>EPAPLT0450B04</td>
      <td>1,5-Naphthalenediamine</td>
      <td>-1.729168</td>
      <td>-2.408105</td>
      <td>36.121572</td>
      <td>35.866971</td>
      <td>74.566382</td>
      <td>74.905851</td>
      <td>73.717711</td>
      <td>-0.686991</td>
      <td>-0.994360</td>
      <td>-0.852628</td>
      <td>16.677642</td>
      <td>12.321542</td>
      <td>14.261381</td>
      <td>-2.068636</td>
      <td>35.994271</td>
      <td>74.396648</td>
      <td>-0.844660</td>
      <td>14.420188</td>
    </tr>
    <tr>
      <th>21</th>
      <td>EPAPLT0458H12</td>
      <td>2,5-Di-tert-butylbenzene-1,4-diol</td>
      <td>-4.275182</td>
      <td>-3.341643</td>
      <td>27.974328</td>
      <td>29.077600</td>
      <td>70.068424</td>
      <td>71.171697</td>
      <td>69.728956</td>
      <td>-4.529101</td>
      <td>0.892544</td>
      <td>-3.328655</td>
      <td>14.073544</td>
      <td>12.903835</td>
      <td>13.397333</td>
      <td>-3.808412</td>
      <td>28.525964</td>
      <td>70.323026</td>
      <td>-2.321738</td>
      <td>13.458238</td>
    </tr>
    <tr>
      <th>0</th>
      <td>EPAPLT0443B12</td>
      <td>Dichlone</td>
      <td>-1.474566</td>
      <td>-1.219965</td>
      <td>67.182942</td>
      <td>64.127725</td>
      <td>55.895613</td>
      <td>56.998886</td>
      <td>59.035697</td>
      <td>10.950335</td>
      <td>11.127926</td>
      <td>11.862196</td>
      <td>22.247848</td>
      <td>22.915522</td>
      <td>25.082472</td>
      <td>-1.347266</td>
      <td>65.655333</td>
      <td>57.310066</td>
      <td>11.313486</td>
      <td>23.415281</td>
    </tr>
    <tr>
      <th>12</th>
      <td>EPAPLT0455G04</td>
      <td>Clofentezine</td>
      <td>-2.238371</td>
      <td>-2.492972</td>
      <td>13.377181</td>
      <td>13.886384</td>
      <td>47.833236</td>
      <td>45.881292</td>
      <td>52.840397</td>
      <td>-3.774340</td>
      <td>-1.771320</td>
      <td>-3.156187</td>
      <td>27.044510</td>
      <td>27.107691</td>
      <td>24.274775</td>
      <td>-2.365671</td>
      <td>13.631783</td>
      <td>48.851642</td>
      <td>-2.900616</td>
      <td>26.142325</td>
    </tr>
    <tr>
      <th>37</th>
      <td>EPAPLT0473E09</td>
      <td>1-Naphthol</td>
      <td>-2.323238</td>
      <td>0.901713</td>
      <td>28.228929</td>
      <td>21.184957</td>
      <td>37.818915</td>
      <td>43.674747</td>
      <td>39.516257</td>
      <td>-5.227512</td>
      <td>-2.802713</td>
      <td>-7.561807</td>
      <td>5.048854</td>
      <td>3.771565</td>
      <td>4.277016</td>
      <td>-0.710762</td>
      <td>24.706943</td>
      <td>40.336640</td>
      <td>-5.197344</td>
      <td>4.365812</td>
    </tr>
    <tr>
      <th>35</th>
      <td>EPAPLT0467F10</td>
      <td>Perfluorooctanoic acid</td>
      <td>-3.002175</td>
      <td>-3.256776</td>
      <td>1.071448</td>
      <td>-1.219965</td>
      <td>14.225853</td>
      <td>9.558161</td>
      <td>45.881292</td>
      <td>-8.787867</td>
      <td>-8.977412</td>
      <td>-8.678581</td>
      <td>125.318868</td>
      <td>123.464410</td>
      <td>123.237298</td>
      <td>-3.129475</td>
      <td>-0.074259</td>
      <td>23.221768</td>
      <td>-8.814620</td>
      <td>124.006859</td>
    </tr>
    <tr>
      <th>38</th>
      <td>EPAPLT0473G07</td>
      <td>5-Chloro-N-(2-chloro-4-nitrophenyl)-2-hydroxyb...</td>
      <td>-3.087042</td>
      <td>-2.747573</td>
      <td>27.719726</td>
      <td>-1.389699</td>
      <td>14.989657</td>
      <td>19.148146</td>
      <td>19.233013</td>
      <td>-7.891375</td>
      <td>-8.125317</td>
      <td>-8.967166</td>
      <td>57.291312</td>
      <td>60.981446</td>
      <td>60.228392</td>
      <td>-2.917308</td>
      <td>13.165014</td>
      <td>17.790272</td>
      <td>-8.327953</td>
      <td>59.500383</td>
    </tr>
    <tr>
      <th>40</th>
      <td>EPAPLT0475F02</td>
      <td>Darbufelone mesylate</td>
      <td>-1.644301</td>
      <td>-2.917308</td>
      <td>-1.135098</td>
      <td>-0.880496</td>
      <td>17.620538</td>
      <td>17.620538</td>
      <td>17.450804</td>
      <td>-8.941552</td>
      <td>-7.597667</td>
      <td>-9.107190</td>
      <td>77.219059</td>
      <td>71.280010</td>
      <td>80.190291</td>
      <td>-2.280804</td>
      <td>-1.007797</td>
      <td>17.563960</td>
      <td>-8.548803</td>
      <td>76.229787</td>
    </tr>
    <tr>
      <th>2</th>
      <td>EPAPLT0444B09</td>
      <td>1-Hydroxypyrene</td>
      <td>-0.625895</td>
      <td>-2.832440</td>
      <td>35.697237</td>
      <td>12.698244</td>
      <td>17.281069</td>
      <td>16.856734</td>
      <td>17.450804</td>
      <td>-9.578488</td>
      <td>-10.007097</td>
      <td>-8.125317</td>
      <td>-3.837521</td>
      <td>-4.071463</td>
      <td>-3.215953</td>
      <td>-1.729168</td>
      <td>24.197740</td>
      <td>17.196202</td>
      <td>-9.236968</td>
      <td>-3.708313</td>
    </tr>
    <tr>
      <th>25</th>
      <td>EPAPLT0462G05</td>
      <td>4-Bromo-2-(4-chlorophenyl)-1-(ethoxymethyl)-5-...</td>
      <td>-2.323238</td>
      <td>-2.747573</td>
      <td>-1.474566</td>
      <td>-2.323238</td>
      <td>26.192118</td>
      <td>10.406832</td>
      <td>14.989657</td>
      <td>-8.617107</td>
      <td>-8.731517</td>
      <td>-9.192570</td>
      <td>17.078929</td>
      <td>24.660694</td>
      <td>21.616034</td>
      <td>-2.535406</td>
      <td>-1.898902</td>
      <td>17.196202</td>
      <td>-8.847064</td>
      <td>21.118552</td>
    </tr>
    <tr>
      <th>1</th>
      <td>EPAPLT0446B07</td>
      <td>7,12-Dimethylbenz(a)anthracene</td>
      <td>-0.371294</td>
      <td>-1.729168</td>
      <td>6.078608</td>
      <td>6.333210</td>
      <td>13.886384</td>
      <td>15.583727</td>
      <td>15.838328</td>
      <td>-5.678319</td>
      <td>-5.715887</td>
      <td>-4.542762</td>
      <td>12.753566</td>
      <td>9.220532</td>
      <td>10.933259</td>
      <td>-1.050231</td>
      <td>6.205909</td>
      <td>15.102813</td>
      <td>-5.312323</td>
      <td>10.969119</td>
    </tr>
    <tr>
      <th>5</th>
      <td>EPAPLT0439F09</td>
      <td>Potassium perfluorooctanesulfonate</td>
      <td>-2.917308</td>
      <td>-3.341643</td>
      <td>0.731979</td>
      <td>-0.625895</td>
      <td>17.026468</td>
      <td>13.716650</td>
      <td>13.971251</td>
      <td>-9.482863</td>
      <td>-8.384873</td>
      <td>-8.135562</td>
      <td>77.582779</td>
      <td>76.889491</td>
      <td>67.062226</td>
      <td>-3.129475</td>
      <td>0.053042</td>
      <td>14.904790</td>
      <td>-8.667766</td>
      <td>73.844832</td>
    </tr>
    <tr>
      <th>29</th>
      <td>EPAPLT0467E10</td>
      <td>Potassium perfluorohexanesulfonate</td>
      <td>-1.814035</td>
      <td>-1.559434</td>
      <td>-1.814035</td>
      <td>-1.729168</td>
      <td>14.056118</td>
      <td>14.395587</td>
      <td>15.498860</td>
      <td>-8.396826</td>
      <td>-8.127024</td>
      <td>-9.021809</td>
      <td>183.894833</td>
      <td>229.648394</td>
      <td>191.702002</td>
      <td>-1.686734</td>
      <td>-1.771601</td>
      <td>14.650188</td>
      <td>-8.515220</td>
      <td>201.748409</td>
    </tr>
    <tr>
      <th>42</th>
      <td>EPAPLT0483C11</td>
      <td>8-Hydroxyquinoline</td>
      <td>-2.917308</td>
      <td>-2.577839</td>
      <td>7.351615</td>
      <td>7.351615</td>
      <td>14.395587</td>
      <td>14.395587</td>
      <td>14.819923</td>
      <td>-4.973079</td>
      <td>-3.945101</td>
      <td>-4.319066</td>
      <td>7.664904</td>
      <td>9.128321</td>
      <td>11.281611</td>
      <td>-2.747573</td>
      <td>7.351615</td>
      <td>14.537032</td>
      <td>-4.412415</td>
      <td>9.358279</td>
    </tr>
    <tr>
      <th>4</th>
      <td>EPAPLT0440A05</td>
      <td>CP-085958</td>
      <td>-2.577839</td>
      <td>-3.087042</td>
      <td>12.104174</td>
      <td>35.527502</td>
      <td>11.000902</td>
      <td>11.170636</td>
      <td>21.439559</td>
      <td>-8.603446</td>
      <td>-6.238414</td>
      <td>-4.353218</td>
      <td>-2.251157</td>
      <td>-4.387370</td>
      <td>-3.559182</td>
      <td>-2.832440</td>
      <td>23.815838</td>
      <td>14.537032</td>
      <td>-6.398359</td>
      <td>-3.399236</td>
    </tr>
    <tr>
      <th>6</th>
      <td>EPAPLT0441F12</td>
      <td>Perfluorooctanesulfonic acid</td>
      <td>-3.681112</td>
      <td>-3.002175</td>
      <td>-1.729168</td>
      <td>15.074524</td>
      <td>12.783111</td>
      <td>14.395587</td>
      <td>15.838328</td>
      <td>-8.743470</td>
      <td>-7.838439</td>
      <td>-7.247608</td>
      <td>196.884582</td>
      <td>194.171198</td>
      <td>197.446384</td>
      <td>-3.341643</td>
      <td>6.672678</td>
      <td>14.339009</td>
      <td>-7.943172</td>
      <td>196.167388</td>
    </tr>
    <tr>
      <th>23</th>
      <td>EPAPLT0458E05</td>
      <td>Napropamide</td>
      <td>-3.596245</td>
      <td>-3.596245</td>
      <td>6.163475</td>
      <td>11.340370</td>
      <td>14.140986</td>
      <td>14.735055</td>
      <td>14.056118</td>
      <td>-6.439911</td>
      <td>-7.541316</td>
      <td>-7.136614</td>
      <td>10.443177</td>
      <td>12.405215</td>
      <td>10.391949</td>
      <td>-3.596245</td>
      <td>8.751923</td>
      <td>14.310720</td>
      <td>-7.039280</td>
      <td>11.080113</td>
    </tr>
    <tr>
      <th>3</th>
      <td>EPAPLT0440C01</td>
      <td>Riboflavin</td>
      <td>2.004986</td>
      <td>0.816846</td>
      <td>28.568398</td>
      <td>7.691084</td>
      <td>14.140986</td>
      <td>13.122580</td>
      <td>15.244258</td>
      <td>-6.434788</td>
      <td>7.068950</td>
      <td>-3.111789</td>
      <td>10.282662</td>
      <td>10.709563</td>
      <td>11.351622</td>
      <td>1.410916</td>
      <td>18.129741</td>
      <td>14.169275</td>
      <td>-0.825876</td>
      <td>10.781283</td>
    </tr>
    <tr>
      <th>16</th>
      <td>EPAPLT0450H08</td>
      <td>Cyclanilide</td>
      <td>-2.662706</td>
      <td>-1.898902</td>
      <td>11.425237</td>
      <td>-0.116692</td>
      <td>13.377181</td>
      <td>13.716650</td>
      <td>14.140986</td>
      <td>-7.568638</td>
      <td>-7.232240</td>
      <td>-7.930650</td>
      <td>140.933206</td>
      <td>151.793572</td>
      <td>153.895633</td>
      <td>-2.280804</td>
      <td>5.654273</td>
      <td>13.744939</td>
      <td>-7.577176</td>
      <td>148.874137</td>
    </tr>
    <tr>
      <th>13</th>
      <td>EPAPLT0450D03</td>
      <td>Perfluorooctanesulfonamide</td>
      <td>-3.596245</td>
      <td>-3.341643</td>
      <td>-1.474566</td>
      <td>-2.153503</td>
      <td>13.377181</td>
      <td>13.546916</td>
      <td>14.056118</td>
      <td>-8.646136</td>
      <td>-8.340475</td>
      <td>-8.673458</td>
      <td>93.965538</td>
      <td>94.496603</td>
      <td>108.292342</td>
      <td>-3.468944</td>
      <td>-1.814035</td>
      <td>13.660072</td>
      <td>-8.553356</td>
      <td>98.918161</td>
    </tr>
    <tr>
      <th>15</th>
      <td>EPAPLT0451B12</td>
      <td>N-Phenyl-1-naphthylamine</td>
      <td>0.562245</td>
      <td>-0.795629</td>
      <td>5.484538</td>
      <td>5.314804</td>
      <td>13.801517</td>
      <td>12.952846</td>
      <td>13.377181</td>
      <td>8.247197</td>
      <td>8.460648</td>
      <td>9.825024</td>
      <td>39.311942</td>
      <td>36.803471</td>
      <td>35.652545</td>
      <td>-0.116692</td>
      <td>5.399671</td>
      <td>13.377181</td>
      <td>8.844290</td>
      <td>37.255986</td>
    </tr>
    <tr>
      <th>17</th>
      <td>EPAPLT0449D03</td>
      <td>meso-Hexestrol</td>
      <td>-2.492972</td>
      <td>-2.747573</td>
      <td>-1.898902</td>
      <td>-2.068636</td>
      <td>13.886384</td>
      <td>12.443643</td>
      <td>12.358776</td>
      <td>-9.569950</td>
      <td>-8.296077</td>
      <td>-9.534091</td>
      <td>51.208824</td>
      <td>54.173226</td>
      <td>57.970939</td>
      <td>-2.620273</td>
      <td>-1.983769</td>
      <td>12.896268</td>
      <td>-9.133373</td>
      <td>54.450996</td>
    </tr>
    <tr>
      <th>45</th>
      <td>EPAPLT0485D02</td>
      <td>Perfluorodecanoic acid</td>
      <td>-3.426510</td>
      <td>-3.171909</td>
      <td>-0.456161</td>
      <td>0.477378</td>
      <td>14.395587</td>
      <td>3.702329</td>
      <td>17.365936</td>
      <td>-8.902277</td>
      <td>-8.717856</td>
      <td>-8.728101</td>
      <td>81.575159</td>
      <td>61.688394</td>
      <td>82.142083</td>
      <td>-3.299210</td>
      <td>0.010608</td>
      <td>11.821284</td>
      <td>-8.782745</td>
      <td>75.135212</td>
    </tr>
    <tr>
      <th>20</th>
      <td>EPAPLT0456G04</td>
      <td>1,4-Benzenediamine</td>
      <td>-3.765979</td>
      <td>-2.747573</td>
      <td>9.133825</td>
      <td>12.443643</td>
      <td>10.491699</td>
      <td>11.594972</td>
      <td>10.746300</td>
      <td>-7.373971</td>
      <td>-2.941029</td>
      <td>-3.487462</td>
      <td>1.020614</td>
      <td>2.217645</td>
      <td>1.167468</td>
      <td>-3.256776</td>
      <td>10.788734</td>
      <td>10.944324</td>
      <td>-4.600821</td>
      <td>1.468576</td>
    </tr>
    <tr>
      <th>18</th>
      <td>EPAPLT0451C01</td>
      <td>8-Hydroxyquinoline sulfate</td>
      <td>-2.068636</td>
      <td>-2.832440</td>
      <td>3.956930</td>
      <td>3.362860</td>
      <td>9.812762</td>
      <td>10.746300</td>
      <td>10.916035</td>
      <td>-6.906087</td>
      <td>-7.397877</td>
      <td>-6.890719</td>
      <td>3.098769</td>
      <td>4.765391</td>
      <td>2.859704</td>
      <td>-2.450538</td>
      <td>3.659895</td>
      <td>10.491699</td>
      <td>-7.064894</td>
      <td>3.574621</td>
    </tr>
    <tr>
      <th>46</th>
      <td>EPAPLT0481H08</td>
      <td>Dinoseb</td>
      <td>-1.983769</td>
      <td>-1.729168</td>
      <td>3.023392</td>
      <td>3.277993</td>
      <td>10.067363</td>
      <td>9.727895</td>
      <td>10.916035</td>
      <td>-6.033501</td>
      <td>-6.731911</td>
      <td>-7.068310</td>
      <td>18.525270</td>
      <td>17.136987</td>
      <td>18.352802</td>
      <td>-1.856468</td>
      <td>3.150692</td>
      <td>10.237098</td>
      <td>-6.611241</td>
      <td>18.005020</td>
    </tr>
    <tr>
      <th>30</th>
      <td>EPAPLT0464E05</td>
      <td>Chlorophacinone</td>
      <td>-1.729168</td>
      <td>-2.323238</td>
      <td>-0.710762</td>
      <td>-1.559434</td>
      <td>8.624622</td>
      <td>9.388426</td>
      <td>9.473293</td>
      <td>-8.321691</td>
      <td>-8.352428</td>
      <td>-8.236311</td>
      <td>21.832900</td>
      <td>20.919332</td>
      <td>21.921696</td>
      <td>-2.026203</td>
      <td>-1.135098</td>
      <td>9.162114</td>
      <td>-8.303477</td>
      <td>21.557976</td>
    </tr>
    <tr>
      <th>34</th>
      <td>EPAPLT0470A10</td>
      <td>Ethylenediaminetetraacetic acid ferric sodium ...</td>
      <td>-2.323238</td>
      <td>-2.153503</td>
      <td>23.306636</td>
      <td>2.174720</td>
      <td>6.333210</td>
      <td>12.698244</td>
      <td>7.860818</td>
      <td>-6.231583</td>
      <td>-4.187580</td>
      <td>-5.396565</td>
      <td>7.401933</td>
      <td>7.473652</td>
      <td>6.394446</td>
      <td>-2.238371</td>
      <td>12.740678</td>
      <td>8.964091</td>
      <td>-5.271909</td>
      <td>7.090010</td>
    </tr>
    <tr>
      <th>33</th>
      <td>EPAPLT0466E09</td>
      <td>5-Amino-2-methylphenol</td>
      <td>-3.087042</td>
      <td>-1.983769</td>
      <td>2.599056</td>
      <td>3.702329</td>
      <td>7.606217</td>
      <td>8.964091</td>
      <td>8.454888</td>
      <td>-8.086042</td>
      <td>-7.672802</td>
      <td>-7.223702</td>
      <td>-1.370033</td>
      <td>-1.364910</td>
      <td>-1.245377</td>
      <td>-2.535406</td>
      <td>3.150692</td>
      <td>8.341732</td>
      <td>-7.660848</td>
      <td>-1.326773</td>
    </tr>
    <tr>
      <th>8</th>
      <td>EPAPLT0446B08</td>
      <td>tert-Butylhydroquinone</td>
      <td>-3.087042</td>
      <td>-3.341643</td>
      <td>61.921180</td>
      <td>2.853657</td>
      <td>8.115419</td>
      <td>8.454888</td>
      <td>8.115419</td>
      <td>-4.626435</td>
      <td>-8.460007</td>
      <td>-8.400241</td>
      <td>-3.326948</td>
      <td>-5.997641</td>
      <td>-2.602923</td>
      <td>-3.214343</td>
      <td>32.387418</td>
      <td>8.228575</td>
      <td>-7.162228</td>
      <td>-3.975837</td>
    </tr>
    <tr>
      <th>32</th>
      <td>EPAPLT0465G08</td>
      <td>1,3-Benzenediamine</td>
      <td>-1.474566</td>
      <td>-0.541028</td>
      <td>3.702329</td>
      <td>3.447727</td>
      <td>8.115419</td>
      <td>7.860818</td>
      <td>8.200286</td>
      <td>-4.851839</td>
      <td>-6.468940</td>
      <td>-6.521876</td>
      <td>1.500451</td>
      <td>1.367258</td>
      <td>1.437270</td>
      <td>-1.007797</td>
      <td>3.575028</td>
      <td>8.058841</td>
      <td>-5.947552</td>
      <td>1.434993</td>
    </tr>
    <tr>
      <th>14</th>
      <td>EPAPLT0448H06</td>
      <td>Nitrofurazone</td>
      <td>-4.699517</td>
      <td>-1.135098</td>
      <td>4.975335</td>
      <td>3.702329</td>
      <td>7.266748</td>
      <td>8.115419</td>
      <td>8.200286</td>
      <td>-5.553664</td>
      <td>-6.286227</td>
      <td>-6.651654</td>
      <td>7.115055</td>
      <td>7.526588</td>
      <td>6.006820</td>
      <td>-2.917308</td>
      <td>4.338832</td>
      <td>7.860818</td>
      <td>-6.163848</td>
      <td>6.882821</td>
    </tr>
    <tr>
      <th>44</th>
      <td>EPAPLT0439F06</td>
      <td>Furazolidone</td>
      <td>-0.031825</td>
      <td>0.562245</td>
      <td>5.314804</td>
      <td>4.126664</td>
      <td>8.030552</td>
      <td>7.521349</td>
      <td>7.860818</td>
      <td>-6.369899</td>
      <td>-6.171817</td>
      <td>-6.774602</td>
      <td>4.814912</td>
      <td>5.964130</td>
      <td>6.242469</td>
      <td>0.265210</td>
      <td>4.720734</td>
      <td>7.804240</td>
      <td>-6.438773</td>
      <td>5.673837</td>
    </tr>
    <tr>
      <th>19</th>
      <td>EPAPLT0447D11</td>
      <td>2,4-Diaminotoluene</td>
      <td>-1.898902</td>
      <td>-1.559434</td>
      <td>6.587811</td>
      <td>3.362860</td>
      <td>6.842412</td>
      <td>7.860818</td>
      <td>7.860818</td>
      <td>-4.891114</td>
      <td>-5.563910</td>
      <td>-6.173525</td>
      <td>3.153412</td>
      <td>4.155776</td>
      <td>2.827260</td>
      <td>-1.729168</td>
      <td>4.975335</td>
      <td>7.521349</td>
      <td>-5.542849</td>
      <td>3.378816</td>
    </tr>
    <tr>
      <th>26</th>
      <td>EPAPLT0463D10</td>
      <td>Nitrofurantoin</td>
      <td>-2.153503</td>
      <td>-2.747573</td>
      <td>2.514189</td>
      <td>2.344454</td>
      <td>5.314804</td>
      <td>5.908874</td>
      <td>5.739140</td>
      <td>-5.821758</td>
      <td>-6.800216</td>
      <td>-6.197431</td>
      <td>5.470632</td>
      <td>6.691569</td>
      <td>6.638633</td>
      <td>-2.450538</td>
      <td>2.429322</td>
      <td>5.654273</td>
      <td>-6.273135</td>
      <td>6.266945</td>
    </tr>
    <tr>
      <th>28</th>
      <td>EPAPLT0466D09</td>
      <td>4-Nitrosodiphenylamine</td>
      <td>2.683923</td>
      <td>2.599056</td>
      <td>12.358776</td>
      <td>-0.201559</td>
      <td>5.739140</td>
      <td>5.739140</td>
      <td>5.229937</td>
      <td>-7.170766</td>
      <td>-8.360966</td>
      <td>-8.099703</td>
      <td>0.329034</td>
      <td>1.435562</td>
      <td>0.885713</td>
      <td>2.641489</td>
      <td>6.078608</td>
      <td>5.569405</td>
      <td>-7.877145</td>
      <td>0.883436</td>
    </tr>
    <tr>
      <th>7</th>
      <td>EPAPLT0439F10</td>
      <td>Phenothiazine</td>
      <td>-4.529783</td>
      <td>-4.529783</td>
      <td>3.362860</td>
      <td>0.222776</td>
      <td>6.163475</td>
      <td>2.344454</td>
      <td>6.842412</td>
      <td>-10.922373</td>
      <td>-9.602395</td>
      <td>-10.329835</td>
      <td>-4.283206</td>
      <td>-3.443065</td>
      <td>-3.292796</td>
      <td>-4.529783</td>
      <td>1.792818</td>
      <td>5.116781</td>
      <td>-10.284868</td>
      <td>-3.673022</td>
    </tr>
    <tr>
      <th>47</th>
      <td>EPAPLT0486F12</td>
      <td>Ethoxyquin</td>
      <td>-1.389699</td>
      <td>-2.068636</td>
      <td>0.477378</td>
      <td>0.392510</td>
      <td>2.768790</td>
      <td>3.702329</td>
      <td>3.872063</td>
      <td>-9.298441</td>
      <td>-9.385529</td>
      <td>-8.383165</td>
      <td>-3.323532</td>
      <td>-3.038362</td>
      <td>-2.985427</td>
      <td>-1.729168</td>
      <td>0.434944</td>
      <td>3.447727</td>
      <td>-9.022379</td>
      <td>-3.115774</td>
    </tr>
    <tr>
      <th>43</th>
      <td>EPAPLT0511B10</td>
      <td>4-Chloro-1,3-diaminobenzene</td>
      <td>-4.105447</td>
      <td>-1.983769</td>
      <td>1.071448</td>
      <td>0.816846</td>
      <td>2.599056</td>
      <td>3.872063</td>
      <td>3.872063</td>
      <td>-7.259561</td>
      <td>-6.678976</td>
      <td>-7.111000</td>
      <td>1.855633</td>
      <td>0.767888</td>
      <td>1.746346</td>
      <td>-3.044608</td>
      <td>0.944147</td>
      <td>3.447727</td>
      <td>-7.016512</td>
      <td>1.456622</td>
    </tr>
    <tr>
      <th>31</th>
      <td>EPAPLT0465B02</td>
      <td>3,5-Dimethyl-1-hexyn-3-ol</td>
      <td>-1.729168</td>
      <td>-1.814035</td>
      <td>-0.710762</td>
      <td>-1.559434</td>
      <td>0.222776</td>
      <td>-0.116692</td>
      <td>-0.116692</td>
      <td>-8.321691</td>
      <td>-8.711025</td>
      <td>-8.705902</td>
      <td>-1.286360</td>
      <td>-1.202687</td>
      <td>-0.446218</td>
      <td>-1.771601</td>
      <td>-1.135098</td>
      <td>-0.003536</td>
      <td>-8.579540</td>
      <td>-0.978422</td>
    </tr>
    <tr>
      <th>27</th>
      <td>EPAPLT0459C02</td>
      <td>alpha-Ionone</td>
      <td>-2.492972</td>
      <td>-4.529783</td>
      <td>-1.219965</td>
      <td>-1.814035</td>
      <td>-0.880496</td>
      <td>0.647112</td>
      <td>-0.371294</td>
      <td>-7.640357</td>
      <td>-7.042695</td>
      <td>-7.339819</td>
      <td>10.555879</td>
      <td>13.738854</td>
      <td>10.309984</td>
      <td>-3.511377</td>
      <td>-1.517000</td>
      <td>-0.201559</td>
      <td>-7.340957</td>
      <td>11.534905</td>
    </tr>
    <tr>
      <th>22</th>
      <td>EPAPLT0463H07</td>
      <td>C10-21 alkanesulfonic acids phenyl esters</td>
      <td>207.468307</td>
      <td>-3.765979</td>
      <td>-2.408105</td>
      <td>-2.153503</td>
      <td>0.222776</td>
      <td>0.222776</td>
      <td>-1.050231</td>
      <td>-8.113364</td>
      <td>-8.641013</td>
      <td>-8.434393</td>
      <td>1.490205</td>
      <td>3.604220</td>
      <td>2.659915</td>
      <td>101.851164</td>
      <td>-2.280804</td>
      <td>-0.201559</td>
      <td>-8.396257</td>
      <td>2.584780</td>
    </tr>
    <tr>
      <th>24</th>
      <td>EPAPLT0458G04</td>
      <td>Imazaquin</td>
      <td>-3.426510</td>
      <td>-3.765979</td>
      <td>-0.880496</td>
      <td>-1.474566</td>
      <td>-0.456161</td>
      <td>-0.965364</td>
      <td>-0.625895</td>
      <td>-6.238414</td>
      <td>-8.815189</td>
      <td>-8.852756</td>
      <td>-1.622758</td>
      <td>-0.196908</td>
      <td>-0.331809</td>
      <td>-3.596245</td>
      <td>-1.177531</td>
      <td>-0.682473</td>
      <td>-7.968786</td>
      <td>-0.717158</td>
    </tr>
    <tr>
      <th>39</th>
      <td>EPAPLT0473A04</td>
      <td>Sulfaquinoxaline</td>
      <td>-3.087042</td>
      <td>-3.087042</td>
      <td>-0.795629</td>
      <td>3.023392</td>
      <td>-1.219965</td>
      <td>-0.116692</td>
      <td>-0.880496</td>
      <td>-8.560756</td>
      <td>-8.519774</td>
      <td>-8.337060</td>
      <td>-1.899390</td>
      <td>-1.202687</td>
      <td>-0.634055</td>
      <td>-3.087042</td>
      <td>1.113881</td>
      <td>-0.739051</td>
      <td>-8.472530</td>
      <td>-1.245377</td>
    </tr>
  </tbody>
</table>
</div>
```
 
 

 
``` python
follow_up_comp_list.to_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ABS_followup_analysis.xlsx")
```
 

 
``` python
follow_up_comp_list_with_mols = follow_up_comp_list.merge(toxcast_library, how='left', on='PREFERRED_NAME').drop_duplicates()
follow_up_comp_list_with_mols = follow_up_comp_list_with_mols.set_index("EPA_SAMPLE_ID")
```
 

 
``` python
follow_up_comp_list_with_mols["SMILES"].fillna("C", inplace=True)
PandasTools.AddMoleculeColumnToFrame(follow_up_comp_list_with_mols, smilesCol='SMILES', molCol='MOL_OBJ', includeFingerprints=True)
```
 

 
``` python
PandasTools.SaveXlsxFromFrame(follow_up_comp_list_with_mols, "C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_followup_raw_analysis.xlsx", molCol='MOL_OBJ')
```
 

 
``` python
```
 
