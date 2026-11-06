"""Independent resources checks; equations preserved from v0.1.11 audit."""

import numpy as np

from .events import close


def vector(value, labels, prefix):
    return (
        np.array([value[f"{prefix}{s + 1}"] for s in labels], float)
        if isinstance(value, dict)
        else np.full(len(labels), value, float)
    )


def resource_shares(host, demand, stock):
    """Reconstruct proportional allocation, including ordered floating-point budget repair."""
    total = np.bincount(host, weights=demand, minlength=len(stock))
    gamma = np.ones(len(stock))
    limited = total > stock
    gamma[limited] = stock[limited] / total[limited]
    x = demand * gamma[host]
    consumed = np.bincount(host, weights=x, minlength=len(stock))

    def ordered_total(values):
        amount = 0.0
        for value in values:
            amount += float(value)
        return amount

    # Sum each over-budget host in consumer order. A common factor preserves equal access;
    # a residual near zero must match this repair, not receive a looser global tolerance.
    for plant in np.flatnonzero(consumed > stock):
        selected = np.flatnonzero(host == plant)
        shares = x[selected].copy()
        used = consumed[plant]
        factor = 1.0
        for _ in range(3):
            if used <= stock[plant]:
                break
            correction = np.nextafter(stock[plant] / used, 0.0)
            factor *= correction
            shares *= correction
            used = ordered_total(shares)
        if used > stock[plant]:
            low, high = 0.0, 1.0
            for _ in range(64):
                middle = (low + high) / 2
                if ordered_total(shares * middle) <= stock[plant]:
                    low = middle
                else:
                    high = middle
            factor *= low
            shares *= low
            used = ordered_total(shares)
        gamma[plant] *= factor
        x[selected] = shares
        consumed[plant] = used
    assert (consumed <= stock).all(), "resource budget"
    close(consumed, np.minimum(total, stock), "proportional resource conservation")
    return total, gamma, x, consumed
