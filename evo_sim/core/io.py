"""Atomic text and table writes, shared by runner and audit."""

import os


def atomic_text(path, text):
    temp = path.with_name(path.name + ".partial")
    temp.write_text(text, encoding="utf-8")
    os.replace(temp, path)


def atomic_csv(frame, path):
    temp = path.with_name(path.name + ".partial")
    frame.to_csv(temp, index=False)
    os.replace(temp, path)
