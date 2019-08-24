/* Arcball, written by Bradley Smith, March 24, 2006
 * arcball.cpp is free to use and modify for any purpose, with no
 * restrictions of copyright or license.
 *
 * See arcball.h for usage details.
 */


#include "arcball.h"
#include <iostream>
#include <fstream>

arcball::arcball()
{
    GLfloat ab_quat[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    GLfloat ab_last[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    GLfloat ab_next[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

    GLdouble ab_glp[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    GLdouble ab_glm[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    int ab_glv[4] = {0,0,640,480};

    //assign values to members
    for(int i=0;i<16;i++)
    {
        m_ab_quat[i]=ab_quat[i];
        m_ab_last[i]=ab_last[i];
        m_ab_next[i]=ab_next[i];
        m_ab_glp[i]=ab_glp[i];
        m_ab_glm[i]=ab_glm[i];
    }

    m_ab_glv[0]=ab_glv[0]; m_ab_glv[1]=ab_glv[1];
    m_ab_glv[2]=ab_glv[2]; m_ab_glv[3]=ab_glv[3];

    // the distance from the origin to the eye
    m_ab_zoom = 1.0;
    m_ab_zoom2 = 1.0;
    // the radius of the arcball
    m_ab_sphere = 1.0;
    m_ab_sphere2 = 1.0;
    // the distance from the origin of the plane that intersects
    // the edge of the visible sphere (tangent to a ray from the eye)
    m_ab_edge = 1.0;
    // whether we are using a sphere or plane
    m_ab_planar = false;
    m_ab_planedist = 0.5;

    m_ab_start = vec(0,0,1);
    m_ab_curr = vec(0,0,1);
    m_ab_eye = vec(0,0,1);
    m_ab_eyedir = vec(0,0,1);
    m_ab_up = vec(0,1,0);
    m_ab_out = vec(1,0,0);
}

void arcball::setzoom(float radius, vec eye, vec up)
{
    m_ab_eye = eye; // store eye vector
    m_ab_zoom2 = m_ab_eye * m_ab_eye;
    m_ab_zoom = sqrt(m_ab_zoom2); // store eye distance
    m_ab_sphere = radius; // sphere radius
    m_ab_sphere2 = m_ab_sphere * m_ab_sphere;
    m_ab_eyedir = m_ab_eye * (1.0 / m_ab_zoom); // distance to eye
    m_ab_edge = m_ab_sphere2 / m_ab_zoom; // plane of visible edge

    if(m_ab_sphere <= 0.0) // trackball mode
    {
        m_ab_planar = true;
        m_ab_up = up;
        m_ab_out = ( m_ab_eyedir ^ m_ab_up );
        m_ab_planedist = (0.0 - m_ab_sphere) * m_ab_zoom;
    }
    else m_ab_planar=false;

    glGetDoublev(GL_PROJECTION_MATRIX,m_ab_glp);
    glGetIntegerv(GL_VIEWPORT,m_ab_glv);
}

void arcball::get_inverse_matrix(float m[16])
{
    m4_inverse(m,m_ab_quat);
}

// affect the arcball's orientation on openGL
void arcball::rotate(void)
{
    glMultMatrixf(m_ab_quat);
}

// convert the quaternion into a rotation matrix
void arcball::quaternion(GLfloat* q, GLfloat x, GLfloat y, GLfloat z, GLfloat w)
{
    GLfloat x2 = x*x;
    GLfloat y2 = y*y;
    GLfloat z2 = z*z;
    GLfloat xy = x*y;
    GLfloat xz = x*z;
    GLfloat yz = y*z;
    GLfloat wx = w*x;
    GLfloat wy = w*y;
    GLfloat wz = w*z;

    q[0] = 1 - 2*y2 - 2*z2;
    q[1] = 2*xy + 2*wz;
    q[2] = 2*xz - 2*wy;

    q[4] = 2*xy - 2*wz;
    q[5] = 1 - 2*x2 - 2*z2;
    q[6] = 2*yz + 2*wx;

    q[8] = 2*xz + 2*wy;
    q[9] = 2*yz - 2*wx;
    q[10]= 1 - 2*x2 - 2*y2;
}

// reset the rotation matrix
void arcball::quatidentity(GLfloat* q)
{
    q[0]=1;  q[1]=0;  q[2]=0;  q[3]=0;
    q[4]=0;  q[5]=1;  q[6]=0;  q[7]=0;
    q[8]=0;  q[9]=0;  q[10]=1; q[11]=0;
    q[12]=0; q[13]=0; q[14]=0; q[15]=1;
}

// copy a rotation matrix
void arcball::quatcopy(GLfloat* dst, GLfloat* src)
{
    dst[0]=src[0]; dst[1]=src[1]; dst[2]=src[2];
    dst[4]=src[4]; dst[5]=src[5]; dst[6]=src[6];
    dst[8]=src[8]; dst[9]=src[9]; dst[10]=src[10];
}

// multiply two rotation matrices
void arcball::quatnext(GLfloat* dest, GLfloat* left, GLfloat* right)
{
    dest[0] = left[0]*right[0] + left[1]*right[4] + left[2] *right[8];
    dest[1] = left[0]*right[1] + left[1]*right[5] + left[2] *right[9];
    dest[2] = left[0]*right[2] + left[1]*right[6] + left[2] *right[10];
    dest[4] = left[4]*right[0] + left[5]*right[4] + left[6] *right[8];
    dest[5] = left[4]*right[1] + left[5]*right[5] + left[6] *right[9];
    dest[6] = left[4]*right[2] + left[5]*right[6] + left[6] *right[10];
    dest[8] = left[8]*right[0] + left[9]*right[4] + left[10]*right[8];
    dest[9] = left[8]*right[1] + left[9]*right[5] + left[10]*right[9];
    dest[10]= left[8]*right[2] + left[9]*right[6] + left[10]*right[10];
}

// find the intersection with the plane through the visible edge
vec arcball::edge_coords(vec m)
{
    // find the intersection of the edge plane and the ray
    float t = (m_ab_edge - m_ab_zoom) / (m_ab_eyedir * m);
    vec a = m_ab_eye + (m*t);
    // find the direction of the eye-axis from that point
    // along the edge plane
    vec c = (m_ab_eyedir * m_ab_edge) - a;

    // find the intersection of the sphere with the ray going from
    // the plane outside the sphere toward the eye-axis.
    float ac = (a*c);
    float c2 = (c*c);
    float q = ( 0.0 - ac - sqrt( ac*ac - c2*((a*a)-m_ab_sphere2) ) ) / c2;

    return (a+(c*q)).unit();
}

// find the intersection with the sphere
vec arcball::sphere_coords(GLdouble mx, GLdouble my)
{
    GLdouble ax,ay,az;

    gluUnProject(mx,my,0,m_ab_glm,m_ab_glp,m_ab_glv,&ax,&ay,&az);
    vec m = vec((float)ax,(float)ay,(float)az) - m_ab_eye;

    // mouse position represents ray: eye + t*m
    // intersecting with a sphere centered at the origin
    GLfloat a = m*m;
    GLfloat b = (m_ab_eye*m);
    GLfloat root = (b*b) - a*(m_ab_zoom2 - m_ab_sphere2);
    if(root <= 0) return edge_coords(m);
    GLfloat t = (0.0 - b - sqrt(root)) / a;
    return (m_ab_eye+(m*t)).unit();
}

// get intersection with plane for "trackball" style rotation
vec arcball::planar_coords(GLdouble mx, GLdouble my)
{
    GLdouble ax,ay,az;

    gluUnProject(mx,my,0,m_ab_glm,m_ab_glp,m_ab_glv,&ax,&ay,&az);
    vec m = vec((float)ax,(float)ay,(float)az) - m_ab_eye;
    // intersect the point with the trackball plane
    GLfloat t = (m_ab_planedist - m_ab_zoom) / (m_ab_eyedir * m);
    vec d = m_ab_eye + m*t;

    return vec(d*m_ab_up,d*m_ab_out,0.0);
}

// reset the arcball
void arcball::reset(void)
{
    quatidentity(m_ab_quat);
    quatidentity(m_ab_last);
}

// begin arcball rotation
void arcball::start(int mx, int my)
{
    // saves a copy of the current rotation for comparison
    quatcopy(m_ab_last,m_ab_quat);
    if(m_ab_planar) m_ab_start = planar_coords((GLdouble)mx,(GLdouble)my);
    else m_ab_start = sphere_coords((GLdouble)mx,(GLdouble)my);
}

// update current arcball rotation
void arcball::move(int mx, int my)
{
    if(m_ab_planar)
    {
        m_ab_curr = planar_coords((GLdouble)mx,(GLdouble)my);
        if(m_ab_curr.equals(m_ab_start)) return;

        // d is motion since the last position
        vec d = m_ab_curr - m_ab_start;

        GLfloat angle = d.length() * 0.5;
        GLfloat cosa = cos( angle );
        GLfloat sina = sin( angle );
        // p is perpendicular to d
        vec p = ((m_ab_out*d.x)-(m_ab_up*d.y)).unit() * sina;

        quaternion(m_ab_next,p.x,p.y,p.z,cosa);
        quatnext(m_ab_quat,m_ab_last,m_ab_next);
        // planar style only ever relates to the last point
        quatcopy(m_ab_last,m_ab_quat);
        m_ab_start = m_ab_curr;
    }
    else
    {
        m_ab_curr = sphere_coords((GLdouble)mx,(GLdouble)my);
        if(m_ab_curr.equals(m_ab_start))
        { // avoid potential rare divide by tiny
            quatcopy(m_ab_quat,m_ab_last);
            return;
        }
        // use a dot product to get the angle between them
        // use a cross product to get the vector to rotate around
        GLfloat cos2a = m_ab_start*m_ab_curr;
        GLfloat sina = sqrt((1.0 - cos2a)*0.5);
        GLfloat cosa = sqrt((1.0 + cos2a)*0.5);
        vec cross = (m_ab_start^m_ab_curr).unit() * sina;
        quaternion(m_ab_next,cross.x,cross.y,cross.z,cosa);
        // update the rotation matrix
        quatnext(m_ab_quat,m_ab_last,m_ab_next);
    }
}

bool m4_inverse(float mr[16],float ma[16])
{
    float  mdet = m4_det(ma);
    float  mtemp[9];

    if (fabs(mdet)<0.0005) return false;

    for (int i=0;i<4;i++)
    for (int j=0;j<4;j++)
    {
        int sign = 1 - ( (i +j) % 2 ) * 2;
        m4_submat( ma, mtemp, i, j );
        mr[i+j*4] = (m3_det(mtemp)*sign)/mdet;
    }

    return true;
}

void m4_submat(float mr[16], float mb[9], int i, int j)
{
    int ti, tj, idst, jdst;

    for ( ti = 0; ti < 4; ti++ )
    {
        if ( ti < i ) idst = ti;
        else if ( ti > i ) idst = ti-1;

        for ( tj = 0; tj < 4; tj++ )
        {
            if ( tj < j ) jdst = tj;
            else if ( tj > j ) jdst = tj-1;

            if ( ti != i && tj != j ) mb[idst*3 + jdst] = mr[ti*4 + tj ];
        }
    }
}

float m4_det(float mr[16])
{
    float  det;
    float  result = 0;
    int    i = 1;
    float  msub3[9];

    for (int n = 0; n < 4; n++, i *= -1 )
    {
        m4_submat( mr, msub3, 0, n );

        det     = m3_det( msub3 );
        result += mr[n] * det * i;
    }

    return result;
}

float m3_det(float mat[9])
{
    float det;

    det = mat[0] * ( mat[4]*mat[8] - mat[7]*mat[5] )
        - mat[1] * ( mat[3]*mat[8] - mat[6]*mat[5] )
        + mat[2] * ( mat[3]*mat[7] - mat[6]*mat[4] );

    return det;
}
