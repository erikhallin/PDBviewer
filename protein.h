#ifndef PROTEIN_H
#define PROTEIN_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <windows.h>
#include <math.h>
#include <map>

#include "arcball.h"
#include "stereo_cam.h"
#include "text_box.h"

#define _piover180 0.0174532925

using namespace std;

enum view_modes{view_arcball,view_locked_roll,view_FPS};

struct atom
{
    float  x,y,z,radius;
    char   type;
    string residue;
};

struct st_block_pos
{
    int block_pos[3];//x,y,z
};

struct st_block_id
{
    int block_id_x,block_id_y,block_id_z;
};

struct st_triangle
{
    int first,second,third;//vertex
    vec first_normal,second_normal,third_normal,
        first_avg_normal,second_avg_normal,third_avg_normal;
};

class protein
{
    public:
        protein();
        ~protein();
        protein(string,HWND,HWND,bool,bool,bool,bool,bool);

        bool m_good;

        operator bool();
        void update_protein(void);
        bool move(bool[],bool[],int[],float);
        void setWindowVal(int,int,int,int);
        void setWindowVal(int,int,int,int,HWND,HWND);

    private:
        HWND m_hWnd_console,m_hWnd_gui; //Windows handles for windows
        text_box Text_box; //for drawing text
        arcball  m_arcball; //for sphere rotation
        //Protein Data
        vector<atom>            m_atoms;
        vector<atom>            m_atoms_other;
        vector< vector<atom> >  m_chains;
        vector<st_block_pos>    m_void_blocks;
        vector<st_block_pos>    m_void_blocks_exp;
        vector<st_block_pos>    m_void_blocks_new;
        vector<st_block_pos>    m_void_blocks_marked;
        vector< vector<vec> >   m_channel_edges_up;
        vector< vector<vec> >   m_channel_edges_down;
        vector< vector<vec> >   m_channel_edges_normals;
        vector< vector<vec> >   m_channel_edges_color;
        vector< vector<float> > m_channel_edges_length;
        vector<vec>             m_mesh_atom_vertices;
        vector<st_triangle>     m_mesh_atom_triangles;
        vector<vec>             m_mesh_channel_vertices;
        vector<st_triangle>     m_mesh_channel_triangles;
        map<int,vector<int> >   m_vertex_triangle_map_atom;
        map<int,vector<int> >   m_vertex_triangle_map_channel;
        string       m_sequence;
        int          m_grid_size[3]; //number of blocks in each dim
        float        m_grid_gap_size;
        bool         ***m_pGrid;
        int          ***m_pGrid_color;
        float        m_solvent_radius;
        float        m_channel_edges_exp_size;
        float        m_channel_edges_step_size;
        int          m_protein_volume;
        float        m_channel_mesh_vol;
        float        m_channel_theta;
        int          m_channel_sidesteps,m_channel_expsteps;
        bool         m_channel_is_leaking,m_channel_defined,m_channel_start_set,m_channel_end_set,m_channel_center_set;
        vec          m_channel_start_pos,m_channel_end_pos,m_channel_center_pos,m_channel_vector;
        //Loading Options
        bool         m_want_mesh_atom,m_want_mesh_channel,m_want_dist_channel,m_want_void_vol,m_want_defined_channel;
        //Controls
        int          m_prev_mouse_pos[2];
        int          m_marked_block[3];
        float        m_axis_roll;
        float        m_key_delay;
        bool         m_marking_mode,m_leakage_waiting_for_response;
        int          m_leakage_option,m_grid_cut_off_left,m_grid_cut_off_right;
        //Viewing
        int          m_view_mode;
        bool         m_anaglyph,m_interlaced,m_switch_eyes,m_mesh_transparency,m_channel_mesh_see_through;
        float        m_eye_separation,m_convergence;
        bool         m_draging,m_moving,m_drag_zooming,m_black_background;
        float        m_rot_hori,m_rot_vert,m_rot_x,m_rot_y,m_rot_z;
        float        m_zoom;
        float        m_center_pos[3],m_min_pos[3],m_max_pos[3];
        int          m_window_height,m_window_length,m_screen_length,m_screen_height;
        float        m_side_shift_view;
        bool         m_show_atoms,m_show_chains,m_show_grid_protein,m_show_grid_void,
                     m_show_channel_dist,m_show_channel_ruler,m_show_mesh_atom,
                     m_show_concave_channel,m_show_other_atoms,m_show_grid_void_exposure,m_show_labels;
        int          m_show_grid_void_exposure_level;
        vec          m_eye,m_eye_dir,m_eye_up,m_eye_right,m_center,m_up,m_axis,m_axis_start,m_axis_end,m_grid_pos,m_light_pos,m_cam_pos;
        //Calculation Flags
        bool         m_tried_to_get_channel_mesh;
        bool         m_calc_done;
        bool         m_grid_updated;
        bool         m_have_mesh_atom;
        bool         m_have_mesh_channel;
        bool         m_channel_edge_set,m_void_vol_set;
        bool         m_have_convex_channel_lengths,m_have_grid_void_color;
        bool         m_have_mesh_normals,m_have_real_mesh_normals;
        bool         m_have_mesh_channel_normals,m_have_real_mesh_channel_normals;
        //Calculation Counters
        int          m_grid_z_progress;
        int          m_axis_roll_progress,m_void_vol_progress;
        int          m_mesh_normal_progress,m_mesh_normal_vertex_progress;

        //Functions
        void         calc_center_pos(void);
        void         calc_grid(void);
        bool         calc_void_volumes_alt(bool);
        bool         calc_grid_void_color(void);
        bool         calc_channel_edges(void);
        bool         calc_convex_channel_lengths(void);
        bool         calc_mesh_normals(void);
        bool         calc_real_mesh_normals(void);
        bool         calc_concave_channel_mesh(void);
        bool         calc_mesh_channel_normals(void);
        bool         calc_real_mesh_channel_normals(void);
        bool         calc_channel_mesh_vol(void);

        void         draw_protein(void);
        void         draw_atoms(void);
        void         draw_alpha_carbons(void);
        void         draw_grid_all(void);
        void         draw_grid_void_volume(void);
        void         draw_grid_void_color(void);
        void         draw_grid_borders(void);
        void         draw_channel_edges(void);
        void         draw_convex_channel(void);
        void         draw_other_atoms(void);
        void         draw_atom_mesh(void);
        void         draw_concave_channel_mesh(void);
        void         draw_labels(void);

        bool         output_diameters_along_axis(void);
        bool         output_diameters_average(void);
        bool         output_diameters_min(void);
        bool         output_radiuses_min(void);

        bool         loadMesh(const char*);
        bool         load_mesh_channel(const char*);
        bool         is_block_inside_planes(st_block_pos);
        void         init_arcball(void);
        void         set_camera_pos(void);
        atom         get_atom_from_string(string);
        bool         new_channel_edges_settings(void);
        bool         collision_water(float,float,float);
};

vec   rotatePointAboutLine(vec p,float theta,vec p1,vec p2);
float volume_of_triangle(vec,vec,vec);
bool  float_to_char(float,char digits[8]);
char  int_2_char(int);
void  get_color(int value,float color[3]);
vec   get_3D_pos(int x, int y);
void  print_controls(void);

#endif // PROTEIN_H
