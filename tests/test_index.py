"""Check shared coupling groups against the alchemical coordinate ordering."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from vaasu.index import write_index_file


def read_index(path):
    groups = {}
    for line in path.read_text().splitlines():
        if line.startswith("["):
            name = line.strip("[] \n")
            groups[name] = []
        elif line.strip() and not line.startswith(";"):
            groups[name].extend(int(value) for value in line.split())
    return groups


class IndexTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        residues = ["GLY", "unk", "POPC", "SOL", "NA", "CL", "B_NA", "B_CL", "COF"]
        records = [
            f"{number:5d}{name:<5}{name:>5}{number:5d}{0.1:8.3f}{0.2:8.3f}{0.3:8.3f}\n"
            for number, name in enumerate(residues, 1)
        ]
        self.coordinates = self.root / "complex_alchemical.gro"
        self.coordinates.write_text("Index fixture\n9\n" + "".join(records) + "2 2 2\n")

    def test_writes_shared_one_based_groups_and_preserves_inputs(self):
        source = self.coordinates.read_bytes()
        working_directory = Path.cwd()
        result = write_index_file(self.root)
        self.assertEqual(result, self.root / "index_alchemical.ndx")
        self.assertEqual(
            read_index(result), {"Solute": [1, 2, 3, 9], "Solvent": [4, 5, 6, 7, 8]}
        )
        self.assertEqual(self.coordinates.read_bytes(), source)
        self.assertEqual(Path.cwd(), working_directory)

    def test_accepts_a_custom_solvent_selection(self):
        result = write_index_file(self.root, solvent_selection="resname SOL")
        self.assertEqual(
            read_index(result), {"Solute": [1, 2, 3, 5, 6, 7, 8, 9], "Solvent": [4]}
        )

    def test_refuses_existing_output_and_dangling_symlinks(self):
        conflict = self.root / "index_alchemical.ndx"
        conflict.write_text("existing result")
        with self.assertRaises(FileExistsError):
            write_index_file(self.root)
        self.assertEqual(conflict.read_text(), "existing result")
        conflict.unlink()
        conflict.symlink_to("missing.ndx")
        with self.assertRaises(FileExistsError):
            write_index_file(self.root)
        self.assertTrue(conflict.is_symlink())

    def test_missing_coordinates_do_not_create_an_index(self):
        self.coordinates.unlink()
        with self.assertRaises(FileNotFoundError):
            write_index_file(self.root)
        self.assertFalse((self.root / "index_alchemical.ndx").exists())


if __name__ == "__main__":
    unittest.main()
