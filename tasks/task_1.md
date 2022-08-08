# Part 1: Data overview

<img width="500px" src="https://habrastorage.org/webt/_k/5u/jh/_k5ujhnchbmrj2d015szelzcyzi.png" />

## Why we consider KPIs?

Procedures for collecting and analyzing of key performance indicators (KPI) are an integral part of modern mobile networks, usually forming part of the  OSS (Operation Support System) functionality. The trend of modern mobile communication standards towards virtualization and softwarization is growing, which is clearly illustrated, for example, by 5G generation networks and studies within the framework of standardization of 6G networks.

Therefore the tasks of developing optimal methods and algorithms that allow, firstly, to process large amounts of data quickly and accurately, and, secondly, automating the planning and optimization of network operation, are one of the most urgent tasks of modern research in the field of communications.


## Main goal of the task

Get the skill of working with tabular data, learn how to make a graphical representation of the data under study in Python.


## Tasks

1. Upload one of the [xlsx-files](https://www.kaggle.com/datasets/vladimirfadeev/lte-technical-kpis) (according to your personal task) into your production environment using the Pandas library tools.

2. Make a primary analysis of the data: what columns, what types of data they contain, the size of the table.

3. Fill in the gaps in the data, if any (some hours and dates of data collection may be missing).

An example of filling for 2017 using the Pandas library:

```python
df = (df.set_index('DT').reindex(pd.date_range(start='1/1/2017', end='31/12/2017', freq='H'))).fillna(method='ffill')
```

4. Visualize the time series

Use `matplotlib` or `seaborn` library (look for `plot()` function).
 
5. Visualize a slice of data over a short period of time.

6. Make a [histogram](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html) of considered data.

7. Build a box-plot of the data.

8. Give a short summary of the analyzed data.


## Working environments options

- Jupyter Notebook (Anaconda)
- Try Jupyter (online)
- Google Colab (online)

## References

1. Fadeev V.A., Zaidullin Sh.V., Fadeeva Z.S., Nadeev A.F. (2021) Monitoring system for LTE-A cellular communication network accessibility indicators. T-Comm, vol. 15, no.3, pр. 4-16. (in Russian)