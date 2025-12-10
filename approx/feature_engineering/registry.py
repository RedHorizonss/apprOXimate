# Allows us to cleanly add feature classes into this dictionary without having to manually change anything.
# Use @registry
FEATURE_REGISTRY = {}

def register_feature(cls):
    """
    Decorator that registers a FeatureModule subclass automatically.
    """
    from .base import FeatureModule
    
    if not issubclass(cls, FeatureModule):
        raise TypeError("Only FeatureModule subclasses can be registered")

    FEATURE_REGISTRY[cls.__name__] = cls
    return cls
