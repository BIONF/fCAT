# fCATpy
Python package for fCAT, a feature-aware completeness assessment tool

python ~/bionf/fCATpy/fcatpy/calcCutoff.py -d test2 -c eukaryota -a eukaryota_busco/weight_dir -b eukaryota_busco/blast_dir --cpus 6 --force

python ~/bionf/fCATpy/fcatpy/searchOrtho.py -d test2 -c eukaryota -a eukaryota_busco/weight_dir -b eukaryota_busco/blast_dir -r "HOMSA@9606@2" -q query/HUMAN@9606@3/HUMAN@9606@3.fa -i 9606 --annoQuery ~/oneseq_data/weight_dir/HUMAN@9606@3.json --cpus 6

python ~/bionf/fCATpy/fcatpy/assessCompleteness.py -d test2 -c eukaryota -o fcatOutput -i HOMSA@9606@1 -m 4
