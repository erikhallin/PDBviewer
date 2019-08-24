#ifndef ARCBALL_H
#define ARCBALL_H

/* Arcball, written by Bradley Smith, March 24, 2006 (modified by me)
 * arcball.h is free to use and modify for any purpose, with no
 * restrictions of copyright or license.
 *
 * Using the arcball:
 *   Call arcball_setzoom after setting up the projection matrix.
 *
 *     The arcball, by default, will act as if a sphere with the given
 *     radius, centred on the origin, can be directly manipulated with
 *     the mouse. Clicking on a point should drag that point to rest under
 *     the current mouse position. eye is the position of the eye relative
 *     to the origin. up is unused.
 *
 *     Alternatively, pass the value: (-radius/|eye|)
 *     This puts the arcball in a mode where the distance the mouse moves
 *     is equivalent to rotation along the axes. This acts much like a
 *     trackball. (It is for this mode that the up vector is required,
 *     which must be a unit vector.)
 *
 *     You should call arcball_setzoom after use of gluLookAt.
 *     gluLookAt(eye.x,eye.y,eye.z, ?,?,?, up.x,up.y,up.z);
 *     The arcball derives its transformation information from the
 *     openGL projection and viewport matrices. (modelview is ignored)
 *
 *     If looking at a point different from the origin, the arcball will still
 *     act as if it centred at (0,0,0). (You can use this to translate
 *     the arcball to some other part of the screen.)
 *
 *   Call arcball_start with a mouse position, and the arcball will
 *     be ready to manipulate. (Call on mouse button down.)
 *   Call arcball_move with a mouse position, and the arcball will
 *     find the rotation necessary to move the start mouse position to
 *     the current mouse position on the sphere. (Call on mouse move.)
 *   Call arcball_rotate after resetting the modelview matrix in your
 *     drawing code. It will call glRotate with its current rotation.
 *   Call arcball_reset if you wish to reset the arcball rotation.
 */

#include "vec.h"
#include <GL\gl.h>
#include <GL\glu.h>

class arcball
{
    public:
        arcball();

        void setzoom(float radius, vec eye, vec up);
        void get_inverse_matrix(float m[16]);
        void rotate(void);
        void rotateb(void);
        void reset(void);
        void start(int mx, int my);
        void move(int mx, int my);

    private:
        GLfloat m_ab_quat[16];
        GLfloat m_ab_last[16];
        GLfloat m_ab_next[16];
        // the distance from the origin to the eye
        GLfloat m_ab_zoom;
        GLfloat m_ab_zoom2;
        // the radius of the arcball
        GLfloat m_ab_sphere;
        GLfloat m_ab_sphere2;
        // the distance from the origin of the plane that intersects
        // the edge of the visible sphere (tangent to a ray from the eye)
        GLfloat m_ab_edge;
        // whether we are using a sphere or plane
        bool m_ab_planar;
        GLfloat m_ab_planedist;
        vec m_ab_start,m_ab_curr,m_ab_eye,m_ab_eyedir,m_ab_up,m_ab_out;
        GLdouble m_ab_glp[16];
        GLdouble m_ab_glm[16];
        int m_ab_glv[4];

        void quaternion(GLfloat* q, GLfloat x, GLfloat y, GLfloat z, GLfloat w);
        void quatidentity(GLfloat* q);
        void quatcopy(GLfloat* dst, GLfloat* src);
        void quatnext(GLfloat* dest, GLfloat* left, GLfloat* right);
        vec  edge_coords(vec m);
        vec  sphere_coords(GLdouble mx, GLdouble my);
        vec  planar_coords(GLdouble mx, GLdouble my);
};

//Matrix inversion
bool  m4_inverse(float mr[16], float ma[16]);
void  m4_submat(float mr[16], float mb[9], int i, int j);
float m4_det(float mr[16]);
float m3_det(float mat[9]);

#endif
