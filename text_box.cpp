#include "text_box.h"

using namespace std;

text_box::text_box()
{
    m_good=true;
    m_dark_background=true;
    //load a font
    load_font(string("data\\Font.bmp"));
}

void text_box::reload(void)
{
    load_font(string("data\\Font.bmp"));
}

void text_box::draw_text(float size,float pos[],string text,vec axisA,vec axisB,float rotA,float rotB)//Locked Roll
{
    glPushMatrix();
    //translate
    glTranslatef(pos[0],pos[1],pos[2]);
    //Locked Roll
    glRotatef(-rotB,axisB.x,axisB.y,axisB.z);
    glRotatef(-rotA,axisA.x,axisA.y,axisA.z);
    glRotatef(-90,0,1,0);
    glRotatef(180,1,0,0);
    //move in front of atom
    glTranslatef(0,0,1);
    //scale
    glScalef(size,size,size);
    draw_text_after_placement(text);
    glPopMatrix();
}

void text_box::draw_text(float size,float pos[],string text,float rota,float rotb)//FPS
{
    glPushMatrix();
    //translate
    glTranslatef(pos[0],pos[1],pos[2]);
    //FPS
    glRotatef(rota+90,0,1,0);
    glRotatef(rotb+90,1,0,0);
    //move in front of atom
    glTranslatef(0,0,1);
    //scale
    glScalef(size,size,size);
    draw_text_after_placement(text);
    glPopMatrix();
}

void text_box::draw_text(float size,float pos[],string text,float matrix[16])//Arcball
{
    glPushMatrix();
    //translate
    glTranslatef(pos[0],pos[1],pos[2]);
    //Arcball
    glMultMatrixf(matrix);
    //move in front of atom
    glTranslatef(0,0,1);
    //scale
    glScalef(size,size,size);
    draw_text_after_placement(text);
    glPopMatrix();
}

void text_box::draw_text_after_placement(string text)
{
    //Draw a transparant square with a letter as texture
    float square[]={-0.38,0.5,0, -0.38,-0.5,0, 0.38,-0.5,0, 0.38,0.5,0};
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glBindTexture(GL_TEXTURE_2D, m_texture_mask_id);
    glVertexPointer(3, GL_FLOAT, 0, square);

    //draw mask
    glDepthMask(GL_FALSE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_DST_COLOR,GL_ZERO);
/*
    //clearing and configuring stencil drawing
    glDrawBuffer(GL_BACK);
    glEnable(GL_STENCIL_TEST);
    glClearStencil(0);
    glClear(GL_STENCIL_BUFFER_BIT);
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE); // colorbuffer is copied to stencil
    glDisable(GL_DEPTH_TEST);
    glStencilFunc(GL_ALWAYS,1,1); // to avoid interaction with stencil content
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
*/

    glPushMatrix();
    for(int letter=0;letter<(int)text.size();letter++)
    {
        //get letter texture offset
        float height_offset=0.72;
        float side_offset;
        char upper_case=text[letter];
        float char_id;
        if(!isupper(upper_case))
        {
            if(!isalpha(upper_case))//number of symbol
            {
                height_offset=0.07;
                if(isdigit(upper_case)) //number
                {
                    char_id=(int)upper_case-48;
                }
                else //symbol
                {
                    char_id=(int)upper_case-32+10;
                }
            }
            else//lower case
            {
                height_offset=0.40;
                upper_case=toupper(upper_case);//make letter upper case for switch
                char_id=(int)upper_case-65;
            }
        }
        else char_id=(int)upper_case-65;//char to int to float

        float letter_tex[]={char_id*0.0385+0.003, height_offset+0.272,
                            char_id*0.0385+0.003, height_offset+0.001,
                            char_id*0.0385+0.040, height_offset+0.001,
                            char_id*0.0385+0.040, height_offset+0.272};
        glTexCoordPointer(2, GL_FLOAT, 0, letter_tex);
        glDrawArrays(GL_QUADS, 0, 4);
        glTranslatef(0.79,0,0);//move to position for next letter
    }
    glPopMatrix();

    //draw text
    glDepthMask(GL_TRUE);
    if(m_dark_background) glBindTexture(GL_TEXTURE_2D, m_texture_id);
    else glBindTexture(GL_TEXTURE_2D, m_texture_dark_id);
    glBlendFunc(GL_ONE,GL_ONE);
    /*glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP); // disabling changes in stencil buffer
    glStencilFunc(GL_EQUAL,1,1); // draws if stencil is 1*/

    //glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    for(int letter=0;letter<(int)text.size();letter++)
    {
        //get letter texture offset
        float height_offset=0.72;
        float side_offset;
        char upper_case=text[letter];
        float char_id;
        if(!isupper(upper_case))
        {
            if(!isalpha(upper_case))//number of symbol
            {
                height_offset=0.07;
                if(isdigit(upper_case)) //number
                {
                    char_id=(int)upper_case-48;
                }
                else //symbol
                {
                    char_id=(int)upper_case-32+10;
                }
            }
            else//lower case
            {
                height_offset=0.40;
                upper_case=toupper(upper_case);//make letter upper case for switch
                char_id=(int)upper_case-65;
            }
        }
        else char_id=(int)upper_case-65;//char to int to float

        float letter_tex[]={char_id*0.0385+0.003, height_offset+0.272,
                            char_id*0.0385+0.003, height_offset+0.001,
                            char_id*0.0385+0.040, height_offset+0.001,
                            char_id*0.0385+0.040, height_offset+0.272};
        glTexCoordPointer(2, GL_FLOAT, 0, letter_tex);
        glDrawArrays(GL_QUADS, 0, 4);
        glTranslatef(0.79,0,0);//move to position for next letter
    }
    glPopMatrix();

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    //glDisable(GL_STENCIL_TEST);
}

void text_box::swap_background_color(void)
{
    m_dark_background=!m_dark_background;
}

bool text_box::load_font(string file)
{
    image* Image;
    Image = loadBMP(file.c_str());
	m_texture_id = loadTexture(Image);

	Image = loadBMP("data\\FontMask.bmp");
	m_texture_mask_id = loadTexture(Image);
	Image = loadBMP("data\\FontDark.bmp");
	m_texture_dark_id = loadTexture(Image);

    return true;
}


//Outside Class

GLuint loadTexture(image* Image) //Makes the image into a texture, and returns the id of the texture
{
	GLuint textureId;
	glGenTextures(1, &textureId); //Make room for our texture
	glBindTexture(GL_TEXTURE_2D, textureId); //Tell OpenGL which texture to edit
	//Map the image to the texture
	glTexImage2D(GL_TEXTURE_2D,                //Always GL_TEXTURE_2D
				 0,                            //0 for now
				 GL_RGB,                       //Format OpenGL uses for image
				 Image->width, Image->height,  //Width and height
				 0,                            //The border of the image
				 GL_RGB, //GL_RGB, because pixels are stored in RGB format
				 GL_UNSIGNED_BYTE, //GL_UNSIGNED_BYTE, because pixels are stored
				                   //as unsigned numbers
				 Image->pixels);               //The actual pixel data
	//Set Texture Settings
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    //GL_NEAREST = pixligt, GL_LINEAR = suddigt

	return textureId; //Returns the id of the texture
}
