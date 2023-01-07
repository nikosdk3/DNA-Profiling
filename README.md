# DNA-Profiling
App that builds DNA strands using a self-implemented vector class and determines who the DNA matches to in the database.

## Usage

### Building
- Clone this repository
- Run `make build`

### Running
- Run `make run`

### Commands
- `load_db`
- `load_dna`
- `process`
- `search`
- `display`

## Example execution

```
Welcome to the DNA Profiling Application.
Enter command or # to exit: load_db small.txt
Loading database...
Enter command or # to exit: display
Database loaded: 
Alice 2 8 3
Bob 4 1 5
Charlie 3 2 5

No DNA loaded.

No DNA has been processed.
Enter command or # to exit: load_dna 1.txt
Loading DNA...
Enter command or # to exit: display
Database loaded: 
Alice 2 8 3
Bob 4 1 5
Charlie 3 2 5

DNA loaded: 
AAGGTAAGTTTAGAATATAAAAGGTGAGTTAAATAGAATAGGTTAAAATTAAAGGAGATCAGATCAGATCAGATCTATCTATCTATCTATCTATCAGAAAAGAGTAAATAGTTAAAGAGTAAGATATTGAATTAATGGAAAATATTGTTGGGGAAAGGAGGGATAGAAGG

No DNA has been processed.
Enter command or # to exit: process
Processing DNA...
Enter command or # to exit: display
Database loaded: 
Alice 2 8 3
Bob 4 1 5
Charlie 3 2 5

DNA loaded: 
AAGGTAAGTTTAGAATATAAAAGGTGAGTTAAATAGAATAGGTTAAAATTAAAGGAGATCAGATCAGATCAGATCTATCTATCTATCTATCTATCAGAAAAGAGTAAATAGTTAAAGAGTAAGATATTGAATTAATGGAAAATATTGTTGGGGAAAGGAGGGATAGAAGG

DNA processed, STR counts: 
AGATC: 4
AATG: 1
TATC: 5

Enter command or # to exit: search
Searching database...
Found in database!  DNA matches: Bob
Enter command or # to exit: #
```
