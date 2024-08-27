```python
import pandas as pd
from glob import glob
import numpy as np
```

## HTS Analysis


```python
filepaths = glob("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\*.csv")
```


```python
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




```python
LOD_avg_cutoff = np.mean(LOD_percentages)
LOD_avg_cutoff
```




    np.float64(1.716948778506576)




```python
LOD_percentages
```




    [np.float64(0.8022782448852165),
     np.float64(1.9631820982061519),
     np.float64(2.3763783935249236),
     np.float64(1.8997665139734345),
     np.float64(5.746864931541303),
     np.float64(0.425204050601918),
     np.float64(2.401031132613852),
     np.float64(0.7711074764966598),
     np.float64(0.9358540780383339),
     np.float64(0.8851304543629532),
     np.float64(1.5925673455582974),
     np.float64(0.8040206222758678)]




```python
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




```python
final_df.sort_values('Rate Average', ascending=False)
```




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




```python
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




```python
final_df = compound_df.join(final_df).sort_values('Rate Average', ascending=False)
```


```python
final_df.plot.scatter(x='R1',y='R2')
```




    <Axes: xlabel='R1', ylabel='R2'>




    
![png](2Y3_HTS_analysis_complete_files/2Y3_HTS_analysis_complete_10_1.png)
    



```python
final_df
```




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




```python
final_LOD_cutoff_df = final_df[final_df["Is_Considered"]]
final_LOD_cutoff_df
```




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




```python
hits_dict = {}
Rate_avg = final_LOD_cutoff_df["Rate Average"]
for i in list(range(100,9,-10))+[5,3,2,1,0]:
    #print(r"Hits >= " + str(i) + r"% of Positive control : ",Rate_avg.where(Rate_avg >= i).count())
    hits_dict[str(i)+r"% and above"] = [Rate_avg.where(Rate_avg >= i).count(), Rate_avg.where(Rate_avg >= i).count()*100/final_df["Rate Average"].count()]
hits_df = pd.DataFrame.from_dict(hits_dict, orient="index", columns=["Number of Hits", "Toxcast Phase II library %"])
hits_df.index.name = "Activity compared to PR"
hits_df
```




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




```python
final_LOD_cutoff_df.to_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_HTS_analysis_above_LOD.xlsx")
```


```python
final_df.count()
```




    EPA_SAMPLE_ID     1920
    PREFERRED_NAME    1920
    R1                1920
    R2                1920
    Is_Considered     1920
    Rate Average      1920
    dtype: int64




```python
final_LOD_cutoff_df.count()
```




    EPA_SAMPLE_ID     110
    PREFERRED_NAME    110
    R1                110
    R2                110
    Is_Considered     110
    Rate Average      110
    dtype: int64




```python
toxcast_library = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\toxcast_library.xlsx")
toxcast_library = toxcast_library[["PREFERRED NAME","SMILES"]]
toxcast_library.rename(columns={"PREFERRED NAME":"PREFERRED_NAME"}, inplace=True)
```


```python
toxcast_library
```




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




```python
final_df_with_mols = final_df.merge(toxcast_library, how='left', on='PREFERRED_NAME').drop_duplicates()
final_df_with_mols = final_df_with_mols.set_index("EPA_SAMPLE_ID")
```


```python
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


```python
final_df_with_mols["SMILES"].fillna("C", inplace=True)
PandasTools.AddMoleculeColumnToFrame(final_df_with_mols, smilesCol='SMILES', molCol='MOL_OBJ', includeFingerprints=True)
```
    


```python
final_df_with_mols_above_lod = final_df_with_mols[final_df_with_mols["Is_Considered"]]
final_df_with_mols_above_lod[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]]
```




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




```python
final_df_with_mols_above_lod[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]].to_html("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\above_LOD_2Y3_screen_data_with_structures.html")
```


```python
final_df_with_mols[["PREFERRED_NAME","R1","R2","Rate Average","MOL_OBJ"]].to_html("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_screen_data_with_structures.html")
```


```python
PandasTools.SaveXlsxFromFrame(final_df_with_mols_above_lod, "C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\above_LOD_2Y3_screen_data_with_structures.xlsx", molCol='MOL_OBJ')

```


```python

```

## HTS Followup


```python
ABS_followup_data = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ABS_followup.xlsx")
ROS_followup_data = pd.read_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ROS_followup.xlsx")
follow_up_picklist = pd.read_csv("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14 Derek_ABS_picklist.csv")
```


```python
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


```python
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


```python
follow_up_comp_list['Avg No Protein'] = follow_up_comp_list[['No Protein R1','No Protein R2']].mean(axis=1)
follow_up_comp_list['Avg Reductase only'] = follow_up_comp_list[['Reductase only R1','Reductase only R2']].mean(axis=1)
follow_up_comp_list['Avg 2Y3-Reductase'] = follow_up_comp_list[['2Y3-Reductase R1','2Y3-Reductase R2','2Y3-Reductase R3']].mean(axis=1)
follow_up_comp_list['Avg Reductase only ROS'] = follow_up_comp_list[['Reductase only ROS R1','Reductase only ROS R2','Reductase only ROS R3']].mean(axis=1)
follow_up_comp_list['Avg 2Y3-Reductase ROS'] = follow_up_comp_list[['2Y3-Reductase ROS R1','2Y3-Reductase ROS R2','2Y3-Reductase ROS R3']].mean(axis=1)
follow_up_comp_list.sort_values("Avg 2Y3-Reductase", ascending=False)
```




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




```python
follow_up_comp_list.to_excel("C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\followup\\2024-02-14_ABS_followup_analysis.xlsx")
```


```python
follow_up_comp_list_with_mols = follow_up_comp_list.merge(toxcast_library, how='left', on='PREFERRED_NAME').drop_duplicates()
follow_up_comp_list_with_mols = follow_up_comp_list_with_mols.set_index("EPA_SAMPLE_ID")
```


```python
follow_up_comp_list_with_mols["SMILES"].fillna("C", inplace=True)
PandasTools.AddMoleculeColumnToFrame(follow_up_comp_list_with_mols, smilesCol='SMILES', molCol='MOL_OBJ', includeFingerprints=True)
```
    

```python
PandasTools.SaveXlsxFromFrame(follow_up_comp_list_with_mols, "C:\\Users\\thisi\\OneDrive\\Desktop\\2Y3_screen_data\\2Y3_followup_raw_analysis.xlsx", molCol='MOL_OBJ')
```


```python

```
