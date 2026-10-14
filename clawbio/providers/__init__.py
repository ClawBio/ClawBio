"""LLM connections shared by skills; task prompts stay with their callers."""

from .chat import GenerationResult, OllamaProvider, OpenAIProvider, ProviderError, create_provider, validate_model_params

__all__ = [
    "GenerationResult",
    "OllamaProvider",
    "OpenAIProvider",
    "ProviderError",
    "create_provider",
    "validate_model_params",
]
