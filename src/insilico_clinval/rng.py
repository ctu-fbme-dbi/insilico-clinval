"""Deterministic pseudo-randomness derived from the suite seed.

All randomness in the framework is a pure function of ``(seed, purpose, key)``:

    derive_int(seed, purpose, key) = int.from_bytes(sha256(f"{seed}|{purpose}|{key}")[:8], "big")

``purpose`` names the decision (e.g. ``"rank"``, ``"haplotype"``) and ``key`` identifies the
object it is made for (e.g. a variant id). Results therefore never depend on iteration order.
"""

import hashlib


def derive_int(seed: int, purpose: str, key: str) -> int:
    """Return a 64-bit unsigned integer derived from the seed, a purpose and a key."""
    digest = hashlib.sha256(f"{seed}|{purpose}|{key}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def derive_choice(seed: int, purpose: str, key: str, n: int) -> int:
    """Return a deterministic index in ``range(n)``."""
    if n < 1:
        raise ValueError("n must be positive")
    return derive_int(seed, purpose, key) % n
