#ifndef TEXT_BOX_H
#define TEXT_BOX_H

#include <iostream>
#include <string>
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <ctype.h>
#include "imageloader.h"
#include "vec.h"

using namespace std;

class text_box
{
    public:
        text_box();

        bool m_good;

        void draw_text(float size,float pos[],string text,vec,vec,float,float);//Locked Roll
        void draw_text(float size,float pos[],string text,float,float);//FPS
        void draw_text(float size,float pos[],string text,float matrix[16]);//arcball
        void draw_text_after_placement(string text);
        void swap_background_color(void);
        void reload(void);

    protected:
    private:
        bool m_dark_background;
        int  m_texture_id;
        int  m_texture_mask_id;
        int  m_texture_dark_id;

        bool load_font(string file);
};

GLuint loadTexture(image* Image); //Makes the image into a texture, and returns the id of the texture

#endif // TEXT_BOX_H
