import csv
import gzip
from itertools import islice

class Utils:
    @staticmethod
    def _to_csv(iterator, file_obj, first_row, chunk_size=5000):
        writer = csv.DictWriter(file_obj, fieldnames=first_row.keys())
        writer.writeheader()
        writer.writerow(first_row)
        while True:
            batch = list(islice(iterator, chunk_size))
            if not batch:
                break
            writer.writerows(batch)
    def to_csv(iterator, filename, chunk_size=5000, buffering=128*1024):
        try:
            first_row = next(iterator)
        except StopIteration:
            return
        
        if filename.endswith('.gz'):
            with gzip.open(filename, 'wt', newline='', buffering=buffering, compresslevel=6, encoding='utf-8') as f:
                Utils._to_csv(iterator, f, first_row, chunk_size)
        else:
            with open(filename, 'w', newline='', buffering=buffering, encoding='utf-8') as f:
                Utils._to_csv(iterator, f, first_row, chunk_size)

