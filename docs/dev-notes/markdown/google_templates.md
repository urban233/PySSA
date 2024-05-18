```python
def process_features(
      self,
      raw_features: Union[tf.train.Example, features.FeatureDict],
      random_seed: int
) -> features.FeatureDict:
    """Processes features to prepare for feeding them into the model.
    
    Args:
      raw_features: The output of the data pipeline either as a dict of NumPy
        arrays or as a tf.train.Example.
      random_seed: The random seed to use when processing the features.
    
    Returns:
      A dict of NumPy feature arrays suitable for feeding into the model.
    """
```
