import os
import chemkin as ck

def check_install():
    return os.path.isdir(ck.ansys_dir)

def test_answer():
    assert check_install()
