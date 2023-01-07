//
//  Project 1 - DNA Profiling
//  Name: Nikos Kotsopulos
//  Program Overview: DNA app that compares DNA strand counts of people loaded from a database against a DNA sequence.
//  Capable of matching a DNA sequence to a person and getting each strand count from the DNA sequence.
//  Creative component: To run creative component type 'getCancerRisk' when prompted and follow steps. Determines a person's cancer risk based on their STR counts.
//
#include <string>
#include <fstream>
#include <sstream>
#include "ourvector.h"

using namespace std;

// Struct that stores info about each person in the database.
struct people
{
    string name;
    ourvector<int> strCount;
};

// Checks to make sure a valid file has been input by the user. Returns 1 is file is valid, 0 if not valid.
int errorHandle(string file)
{
    ifstream inFS;

    inFS.open(file);
    if (!inFS.is_open())
    {
        cout << "Error: unable to open '" << file << "'" << endl;
        return 0;
    }
    return 1;
}

// Function to read in file and checks if valid db was loaded.
int loadDB(string &file)
{
    ifstream inFS(file);
    string data, line;
    while (getline(inFS, data, '\n'))
    {
        line += data + '\n'; // pushing everything into string
    }
    file = line;
    inFS.close();
    return 1;
}

// pushes back the DNA strands into 2d char vector strseq to access entire strand and individual characters for later.
void storeSeqTypes(string &file, ourvector<ourvector<char>> &strseq)
{
    ourvector<char> tempDNA;
    tempDNA.clear(); // clear in case new file is read in.
    strseq.clear();  // clear in case new file is read in.
    long unsigned int index = file.find(",") + 1;

    for (; index < (file.substr(index, file.find('\n'))).size(); index++)
    {
        if (file[index] == ',') // once , is found store all chars into strseq
        {
            strseq.push_back(tempDNA);
            tempDNA.clear();
        }
        else
        {
            tempDNA.push_back(file[index]); // store characters
        }
    }
    strseq.push_back(tempDNA); // push back last str
    tempDNA.clear();
    file.erase(0, file.find('\n')); // erase so that function getSTRCounts can start reading the str counts.
}

// Function that reads the STR counts of each individual and stores them into strcount
void getSTRCounts(string &file, ourvector<int> &strcount)
{
    strcount.clear();
    int convertInt = 0;
    string convertString = "";
    int strandSize = file.substr(0, file.find('\n')).size();
    for (int i = 0; i < strandSize + 1; i++)
    {
        if (file.at(i) == ',' || file.at(i) == '\n') // don't want to push_back , or \n
        {
            convertInt = stoi(convertString); // conversion from string to int to store in ourvector<int>
            convertString = "";
            strcount.push_back(convertInt);
        }
        else
        {
            convertString.push_back(file.at(i)); // in case of double digit number store as string first and then once , or \n is found convert to int and store.
        }
    }
    file.erase(0, strandSize);
}

// stores data from previous two functions into struct profiles
void storeData(string &file, ourvector<int> &strcount, ourvector<people> &profiles)
{
    people newProfile; // declare instance of struct.
    string newName;
    profiles.clear();    // clear in case another database is loaded in.
    while (file != "\n") // looping until \n which means new person
    {
        int strpos = file.find(",");
        newName = file.substr(0, strpos);
        newProfile.name = newName; // push name into newPeople
        file.erase(0, strpos + 1); // erase before calling getSTRCounts since I don't want the name to be read while getting counts
        getSTRCounts(file, strcount);
        newProfile.strCount = strcount; // push strCounts into newPeople
        profiles.push_back(newProfile); // pushes entire newPeople into profiles struct.
    }
    file.clear(); // clear to prepare for potential new file read
}

//Displays database from info stores in the vector of structs.
void displayDatabase(ourvector<people> &profiles)
{
    cout << "Database loaded: ";
    for (int i = 0; i < profiles.size(); i++)
    {
        int strSize = profiles.at(i).strCount.size();
        cout << profiles.at(i).name << " ";
        for (int j = 0; j < strSize; j++)
        {
            if (j < strSize - 1)
            {
                cout << profiles.at(i).strCount.at(j) << " ";
            }
            if (j == strSize - 1)
            {
                cout << profiles.at(i).strCount.at(j);
            }
        }
    }
    cout << endl;
}

//loading DNA from DNA file and pushing back each char into vec.
void loadDNA(string file, ourvector<char> &vec)
{
    vec.clear();//clear to make sure another DNA file does not get loaded
    ifstream inFS(file);
    string data;
    inFS >> data;
    for (auto e : data)
    {
        vec.push_back(e);//pushing back each char.
    }
}

//Display each char from vec
void displayLoadedDNA(ourvector<char> &vec)
{
    cout << endl
         << "DNA Loaded:" << endl;
    for (int i = 0; i < vec.size(); i++)
    {
        cout << vec.at(i);
    }
    cout << endl;
}

// Helper function to process DNA. gets loop values from processDNA and stores them into DNAindex and strIndex. Uses 
// these values to compare each str in strseq to the DNA sequence loaded. If chars don't match the loop break and goes back to the
// processDNA function. If str is found in DNA, it is used in processDNA function to see if the sizes match.
void processHelper(ourvector<ourvector<char>> &strseq, ourvector<char> &dnaseq, int &DNAindex, int &strCount, int strIndex)
{
    for (int i = 0; i < strseq.at(strIndex).size(); i++)
    {
        if (strseq.at(strIndex).at(i) != dnaseq.at(DNAindex))
        {
            strCount = 0;
            break;
        }
        else
        {
            DNAindex++;
            if (DNAindex == dnaseq.size())
            {
                break;
            }
            strCount++;
        }
    }
}

//processes DNA. goes through every strand in the database and checks it against every index in the DNA sequence. longest strand sequence
//gets saved in maxcount and this is what gets processed from the vector strCount.
void processDNA(ourvector<ourvector<char>> &strseq, ourvector<people> &allprof, ourvector<char> &dnaseq, ourvector<int> &strCount)
{
    cout << "Processing DNA..." << endl;
    int strcount = 0, STRrepeats = 0, maxcount = 0, strIndex = 0;
    strCount.clear();
    for (int i = 0; i < dnaseq.size(); i++)//loop for how many different types of strs
    {
        int DNAIndex = i;
        processHelper(strseq, dnaseq, DNAIndex, strcount, strIndex);
        if(strcount%strseq.at(strIndex).size()==0 && strcount!=0){
            STRrepeats++;
            i = DNAIndex - 1;
        } else{//once its not consecutive anymore save the str count to maxcount 
            if(maxcount<STRrepeats){
                maxcount = STRrepeats;
            }
            STRrepeats=0;//resetting counter
        }
        if (i == (dnaseq.size() - 1) && strIndex < strseq.size())//reaches end and resetting
        {
            strCount.push_back(maxcount);
            strcount = 0;
            i=-1;
            strIndex++;
            maxcount = 0;
            STRrepeats = 0;

            if (strIndex == strseq.size())
            {
                break;
            }
        }
    }
}

//displays strs and strCount matches in a formatted way.
void displayProcessedDNA(ourvector<ourvector<char>> &strseq, ourvector<int> &strCount)
{
    cout << endl
         << "DNA processed, STR counts: " << endl;
    for (int i = 0; i < strseq.size(); i++)
    {
        for (int j = 0; j < strseq.at(i).size(); j++)
        {
            cout << strseq.at(i).at(j);
        }
        cout << ": " << strCount.at(i) << endl;
    }
    cout << endl;
}

//uses the processed DNA to find a match from the str counts in the struct.
void searchDB(ourvector<people> &foundPerson, ourvector<int> &foundSTR)
{
    cout << "Searching Database..." << endl;
    int strIndex = 0;
    for (int i = 0; i < foundPerson.size(); i++)
    {
        for (int j = 0; j < foundSTR.size(); j++)
        {
            if (foundPerson.at(i).strCount.at(j) == foundSTR.at(j))//check if strand count in struct is equal to what was processed
            {
                strIndex++;
            } 
            else//make sure all strands for a specific person match
            {
                strIndex = 0;
                break;
            }
        }
        if (strIndex == foundSTR.size())
        {
            foundPerson.at(i).name.erase(0, 1);//to get rid of \n char
            cout << "Found in database!  DNA matches: " << foundPerson.at(i).name << endl;
            foundPerson.at(i).name.insert(0,"\n");
            return;
        }
    }
    cout << "Not found in database." << endl;
}


//Displays all data depending on if flags are true or false. Flags are true or false depending on if certain commands have been ran before.
void displayData(ourvector<ourvector<char>> &strseq, ourvector<people> &profs, ourvector<char> &dnaseq, ourvector<int> &strCount, bool flag1, bool flag2, bool flag3)
{
    if (flag1)
    {
        displayDatabase(profs);
    }
    else
    {
        cout << "No database loaded." << endl;
    }
    if (flag2)
    {
        displayLoadedDNA(dnaseq);
    }
    else
    {
        cout << endl
             << "No DNA loaded." << endl;
    }
    if (flag3)
    {
        displayProcessedDNA(strseq, strCount);
    }
    else
    {
        cout << endl
             << "No DNA has been processed." << endl;
    }
}

//helper function for getCancerRisk. Uses people as a parameter and gets the max str size in the current loaded database to use
//as modulo for rand.
int getMaxStrCount(ourvector<people> &people)
{
    int max = 0;
    for (int i = 0; i < people.size(); i++)
    {
        for (int j = 0; j < people.at(i).strCount.size(); j++)
        {
            if (people.at(i).strCount.at(j) > max)
            {
                max = people.at(i).strCount.at(j);
            }
        }
    }
    return max;
}

//Checks if the individual the user inputs is loaded in the people struct/in the database.
bool getValidIndiviudal(ourvector<people> &ppl, string name, int &index)
{
    for (index = 0; index < ppl.size(); index++)
    {
        if (name == ppl.at(index).name)
        {
            return true;
        }
    }
    return false;
}

// Gets cancer risk for a given indivudal
// If a specific STR count matches with the badSTR count which is randomly generated, the cancer risk for the individual increases.
// Sorta like the lottery except you get unluckier the more the STR and badSTR counts match up.
void getCancerRisk(ourvector<people> &ppl, ourvector<int> &strcount)
{
    string name;
    ourvector<int> badStr;//stores the bad str in here.
    double cancerRisk = 0;
    int index = 0;
    cout << "Enter individual's name: ";
    cin >> name;
    string checkinp = "\n" + name;//adding \n before the name since this is how it was stored in my struct.
    if (getValidIndiviudal(ppl, checkinp, index))//first check if input is valid
    {
        for (int j = 0; j < ppl.at(index).strCount.size(); j++)//
        {
            int badStrCount = rand() % getMaxStrCount(ppl);//generates a random strcount to use as bad strand.
            badStr.push_back(badStrCount);
            if (ppl.at(index).strCount.at(j) == badStr.at(j))
                cancerRisk += (1 / (double)(ppl.at(index).strCount.size())) * 100;//increase cancer risk by 1/(# of different types of strs)*100 each time 
        }
        cout << name << "'s cancer risk is " << cancerRisk << "%" << endl;
        cout << "The bad strand count was ";
        for (int i = 0; i < badStr.size(); i++)
        {
            cout << badStr.at(i) << " ";
        }
        cout << endl;
    }
    else
        cout << "Invalid individual." << endl;
}

//handles error checking before running process, search, or getCancerRisk. returns true or false depending on if database or 
// dna has been loaded or processed.
bool handleErrors(string input, bool flag1, bool flag2, bool flag3)
{
    if ((input == "process" || input == "search" || input == "getCancerRisk")&&!flag1){
        cout << "No database loaded."<<endl;
        return 0;
    }
    if ((input == "process" || input == "search")&&!flag2){
        cout << "No DNA loaded." << endl;
        return 0;
    }
    if(input =="search"&&!flag3){
        cout << "No DNA processed." << endl;
        return 0;
    }
        return 1;
}

//The "driver" function of the program. Takes in all of the user input and runs the appropriate 
//functions to make the code work correctly.
void selectChoice(string selectChoice, string &dbFile, string &dnaFile, ourvector<char> &vec, ourvector<int> &strOccurences, ourvector<int> &strPerPerson, ourvector<people> &profiles, ourvector<ourvector<char>> &seq, bool &flag1, bool &flag2, bool &flag3)
{
    if (selectChoice == "load_db"){//loading database
        cin >> dbFile;//database file
        cout << "Loading database..." << endl;
        if (errorHandle(dbFile)){//checks if database file is valid
            flag1 = 1;//if valid set respective flag to true.
            //runs appropriate functions
            loadDB(dbFile);
            storeSeqTypes(dbFile, seq);
            storeData(dbFile, strPerPerson, profiles);
        } else
            flag1 = 0;
    }
    else if (selectChoice == "load_dna")//loading dna
    {
        cin >> dnaFile; //dna file
        cout << "Loading DNA..." << endl;
        if (errorHandle(dnaFile)){ //checks if dna file is valid
            flag2 = 1;
            loadDNA(dnaFile, vec);
        } else
            flag2 = 0;
    }
    else if (selectChoice == "process"){
        if (handleErrors(selectChoice, flag1, flag2, flag3))//checks if previous commands were ran successfully.
        {
            flag3 = 1;
            processDNA(seq, profiles, vec, strOccurences);
        }
    }
    else if (selectChoice == "search"){
        if (handleErrors(selectChoice, flag1, flag2, flag3))
            searchDB(profiles, strOccurences);
    }
    else if (selectChoice == "getCancerRisk"){
        if (handleErrors(selectChoice, flag1, flag2, flag3))
            getCancerRisk(profiles, strOccurences);
    }
    else if (selectChoice == "display")
        displayData(seq, profiles, vec, strOccurences, flag1, flag2, flag3);
}

int main()
{
    srand(time(NULL)); //used for creative component.
    string userChoice, databaseFile, DNAfile;
    ourvector<char> dnaVec;//vector for dna sequence.
    ourvector<int> STRinDNAcount, strsPerPerson;//vector for matches in the dna, vector for how many unique strs in the database
    ourvector<people> dnaProfiles; //vector of structs used to store name and str count info.
    ourvector<ourvector<char>> sequences; //2d vector to use for process. Stores each str from the database and can access each char at each str.
    bool dbflag = 0, DNAflag = 0, processDNAFlag = 0;
    cout << "Welcome to the DNA Profiling Application." << endl;
    while (userChoice != "#")
    {
        cout << "Enter command or # to exit: ";
        cin >> userChoice;
        selectChoice(userChoice, databaseFile, DNAfile, dnaVec, STRinDNAcount, strsPerPerson, dnaProfiles, sequences, dbflag, DNAflag, processDNAFlag);
    }
    return 0;
}