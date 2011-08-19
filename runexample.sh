echo -e 'Ih Lennard-Jones clusters'
echo -e 'N\tQ_6\tQ_10\tW_6\tW_10'
./example.out 13 13 1.25 6 10
./example.out 55 55 1.25 6 10
./example.out 147 147 1.25 6 10
echo -e '\n'
echo -e 'Oh Morse cluster'
echo -e 'N\tQ_2\tQ_4\tQ_6\tQ_8\tQ_10\tW_2\tW_4\tW_6\tW_8\tW_10'
./example.out 38DD 38 1. 2 4 6 8 10
echo -e '\n'
echo -e 'Td Lennard-Jones clusters'
echo -e 'N\tQ_2\tQ_4\tQ_6\tQ_8\tQ_10\tW_2\tW_4\tW_6\tW_8\tW_10'
./example.out 4 4 1.25 2 4 6 8 10
./example.out 26 26 1.25 2 4 6 8 10
./example.out 98 98 1.25 2 4 6 8 10
echo -e '\n'
echo -e 'D5h Morse clusters'
echo -e 'N\tQ_2\tQ_4\tQ_6\tQ_8\tQ_10\tW_2\tW_4\tW_6\tW_8\tW_10'
./example.out 7AA 7 1.1 2 4 6 8 10
./example.out 13BB 13 1.2 2 4 6 8 10
./example.out 19AA 19 1.2 2 4 6 8 10
./example.out 75CC 75 1.2 2 4 6 8 10
echo -e '\n'
echo -e '39 atom Morse clusters'
echo -e 'N\tQ_2\tQ_4\tQ_6\tQ_8\tQ_10\tW_2\tW_4\tW_6\tW_8\tW_10'
./example.out 39AA 39 1.2 2 4 6 8 10
./example.out 39BB 39 1.2 2 4 6 8 10
./example.out 39CC 39 1.2 2 4 6 8 10
echo -e '\n'