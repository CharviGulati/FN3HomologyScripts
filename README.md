# Python Scripts

## The following commands show how to run specific scripts and how to save their logs to a log file

To run PSIBlastScript.py:

```python3 PSIBlastScript.py 2>&1 | tee LogFiles/psiblast_run.log```

To run FN3_MethodFiltering.py

```python3 FN3_MethodFiltering.py 2>&1 | tee LogFiles/filtering.log```

To run AdnectinHomologPDBScript.py:

```python3 AdnectinHomologPDBScript.py 2>&1 | tee LogFiles/psiblast_filter_and_alignment.log```

To run FN3InterProData.py:

```python3 FN3InterProData.py 2>&1 | tee LogFiles/interpro.log```

To run FN3_merge_align.py:

```python3 FN3_merge_align.py 2>&1 | tee LogFiles/test.log```
