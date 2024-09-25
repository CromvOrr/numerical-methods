CC := g++
FLAGS := -o

SORT := "sort.cpp"
CSORT := "sort.exe"

SRC_A := "num-01-a.cpp"
CSRC_A := "num-01-a.exe"
SRC_B := "num-01-b.cpp"
CSRC_B := "num-01-b.exe"

FLOAT_FILE_A := "data_float_a.dat"
DOUBLE_FILE_A := "data_double_a.dat"
FLOAT_SORTED_FILE_A := "data_float_sorted_a.dat"
DOUBLE_SORTED_FILE_A := "data_double_sorted_a.dat"

FLOAT_FILE_B := "data_float_b.dat"
DOUBLE_FILE_B := "data_double_b.dat"
FLOAT_SORTED_FILE_B := "data_float_sorted_b.dat"
DOUBLE_SORTED_FILE_B := "data_double_sorted_b.dat"

run:
	$(CC) $(SRC_A) $(FLAGS) $(CSRC_A)
	.\$(CSRC_A)
	$(CC) $(SRC_B) $(FLAGS) $(CSRC_B)
	.\$(CSRC_B)
	$(CC) $(SORT) $(FLAGS) $(CSORT)
	.\$(CSORT) $(FLOAT_FILE_A) $(FLOAT_SORTED_FILE_A)
	.\$(CSORT) $(DOUBLE_FILE_A) $(DOUBLE_SORTED_FILE_A)
	.\$(CSORT) $(FLOAT_FILE_B) $(FLOAT_SORTED_FILE_B)
	.\$(CSORT) $(DOUBLE_FILE_B) $(DOUBLE_SORTED_FILE_B)

clean:
	del *.dat *.exe