# gaussprot
Utility for generation Gaussian kernel models of proteins 

```shell
pip install gaussprot
```


```python
from gaussprot import GaussProt

schema = {'A': 1, 'G': '4', 'P': '3', 'R': 1}

gp = GaussProt(schema, standardize_schema=True, model_type='discrete')

models = gp.generate_models(['LPEREDF', 'AQWDCFV'])
```

