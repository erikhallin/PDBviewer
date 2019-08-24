#include "menu.h"

menu::menu()
{
    m_got_filename=false;
    m_want_mesh_atom=m_want_mesh_channel=m_want_dist_channel=m_want_void_vol=true;
    m_want_defined_channel=false;
}

int menu::menu_input(bool keys[])
{
    if (keys[13])//done
    {
        m_got_filename=true;
        //test if file exists
        ifstream file(m_filename.c_str());
        if(file) return 1; //file found
        else return 2; //file not found
    }

    //get list of pressed keys
    int pressed_key=0;
    for(int index=20;index<122;index++)
    {
        if(keys[index])
        {
            //check if keys are valid
            if(index>=21 && index<=47) continue;
            if(index>=58 && index<=64) continue;
            if(index>=91 && index<=94) continue;
            if(index==96) continue;
            pressed_key=index;//store
            break; //dont check more than one key
        }
    }
    //store new char
    m_filename.append(1,int_to_char(pressed_key));
    return 0;
}

string menu::get_filename(void)
{
    return m_filename;
}

bool menu::set_filename(string sInput)
{
    //remove loading options
    istringstream iss(sInput);//for deviding path from settings
    string word;
    iss>>word;//get path
    m_filename=word;
    iss>>word;//get if want atom mesh
    if(word=="0") m_want_mesh_atom=false; else m_want_mesh_atom=true;
    iss>>word;//get if want channel mesh
    if(word=="0") m_want_mesh_channel=false; else m_want_mesh_channel=true;
    iss>>word;//get if want channel volume
    if(word=="0") m_want_void_vol=false; else m_want_void_vol=true;
    iss>>word;//get if want channel dist
    if(word=="0") m_want_dist_channel=false; else m_want_dist_channel=true;
    iss>>word;//get if want to define channel
    if(word=="1") m_want_defined_channel=true; else m_want_defined_channel=false;
    //channel mesh require void vol
    if(m_want_mesh_channel) m_want_void_vol=true;

    ifstream test(m_filename.c_str());
    if(test==0)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"\n***********************\nCould not find that file, try again.\n***********************\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    return true;
}

char int_to_char(int i) //input key to char conversion
{
    switch(i)
    {
        case 20: return ' ';
        case 48: return '0';
        case 49: return '1';
        case 50: return '2';
        case 51: return '3';
        case 52: return '4';
        case 53: return '5';
        case 54: return '6';
        case 55: return '7';
        case 56: return '8';
        case 57: return '9';
        case 65: return 'A';
        case 66: return 'B';
        case 67: return 'C';
        case 68: return 'D';
        case 69: return 'E';
        case 70: return 'F';
        case 71: return 'G';
        case 72: return 'H';
        case 73: return 'I';
        case 74: return 'J';
        case 75: return 'K';
        case 76: return 'L';
        case 77: return 'M';
        case 78: return 'N';
        case 79: return 'O';
        case 80: return 'P';
        case 81: return 'Q';
        case 82: return 'R';
        case 83: return 'S';
        case 84: return 'T';
        case 85: return 'U';
        case 86: return 'V';
        case 87: return 'W';
        case 88: return 'X';
        case 89: return 'Y';
        case 90: return 'Z';
        case 95: return '_';
        case 97: return 'a';
        case 98: return 'b';
        case 99: return 'c';
        case 100: return 'd';
        case 101: return 'e';
        case 102: return 'f';
        case 103: return 'g';
        case 104: return 'h';
        case 105: return 'i';
        case 106: return 'j';
        case 107: return 'k';
        case 108: return 'l';
        case 109: return 'm';
        case 110: return 'n';
        case 111: return 'o';
        case 112: return 'p';
        case 113: return 'q';
        case 114: return 'r';
        case 115: return 's';
        case 116: return 't';
        case 117: return 'u';
        case 118: return 'v';
        case 119: return 'w';
        case 120: return 'x';
        case 121: return 'y';
        case 122: return 'z';
    }
    return 'x';//error
}

void menu::draw_menu(void)
{
    //not used
}

void menu::reset(void)
{
    m_filename.erase();
    m_got_filename=false;
}
