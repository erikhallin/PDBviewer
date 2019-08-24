#ifndef MENU_H
#define MENU_H

#include <windows.h>
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>


using namespace std;

char int_to_char(int);

class menu
{
    public:
        menu();

        bool    m_want_mesh_atom,m_want_mesh_channel,m_want_dist_channel,m_want_void_vol,m_want_defined_channel;

        int     menu_input(bool[]);
        string  get_filename(void);
        bool    set_filename(string);
        void    draw_menu(void);
        void    reset(void);
    private:
        string  m_filename;
        bool    m_got_filename;
};

#endif // MENU_H
