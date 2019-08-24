#include "protein.h"

/*TO-DO

-separat process som laddar mesh?

-är är pol 3d fel öga?

-ruler pekar inte upp vid vissa kanal vec

*/



using namespace std;

protein::protein()
{
    m_good=false;
}

protein::~protein()
{
    //clean up
    delete[] m_pGrid;
    delete[] m_pGrid_color;
    m_pGrid=0;
    m_pGrid_color=0;
    m_atoms.clear();
    m_atoms_other.clear();
    m_chains.clear();
    m_void_blocks.clear();
    m_void_blocks_exp.clear();
    m_void_blocks_new.clear();
    m_channel_edges_up.clear();
    m_channel_edges_down.clear();
    m_channel_edges_normals.clear();
    m_channel_edges_color.clear();
    m_channel_edges_length.clear();
    m_mesh_atom_vertices.clear();
    m_mesh_atom_triangles.clear();
    m_mesh_channel_vertices.clear();
    m_mesh_channel_triangles.clear();
    m_vertex_triangle_map_atom.clear();
    m_vertex_triangle_map_channel.clear();
}

protein::protein(string file_name,HWND hWnd_gui,HWND hWnd_console,bool want_mesh_atom,bool want_mesh_channel,bool want_void_vol,
                 bool want_dist_channel,bool want_defined_channel)
{
    //Save windows handles
    m_hWnd_gui=hWnd_gui;
    m_hWnd_console=hWnd_console;
    //Loading options
    m_want_mesh_atom=want_mesh_atom;
    m_want_mesh_channel=want_mesh_channel;
    m_want_dist_channel=want_dist_channel;
    m_want_void_vol=want_void_vol;
    m_want_defined_channel=want_defined_channel;
    m_channel_defined=false;
    if(want_defined_channel) want_dist_channel=false;//not compatible
    if(!want_void_vol) m_want_mesh_channel=false;//channel mesh requires void vol

    ifstream file(file_name.c_str());
    if (file==0)
    {
        m_good=false;//Bad file
    }
    else //File OK
    {
        //Read positions and type of atoms
        enum chains{chainA=0,chainB,chainC,chainD,chainE,chainF,
                    chainG,chainH,chainI,chainJ,chainK,chainL};
        char curr_chain='w'; int curr_chain_ID=-1;
        string line;
        while (getline(file,line)) //read PDB file
        {
            if (line[0]=='A'&&line[1]=='T'&&line[2]=='O'&&line[3]=='M')
            {
                if(line[13]!='C'&&line[13]!='O'&&line[13]!='N'&&line[13]!='S'&&line[13]!='H')//other atom
                {
                    atom temp_atom=get_atom_from_string(line);
                    m_atoms_other.push_back(temp_atom);
                }
                else//protein atom
                {
                    atom temp_atom=get_atom_from_string(line);
                    m_atoms.push_back(temp_atom);
                    if(line[13]=='C' && line[14]=='A') //atom is alpha carbon
                    {
                        if(curr_chain!=line[21])
                        {
                            curr_chain=line[21];

                            m_chains.push_back(vector<atom>());
                            curr_chain_ID++;
                        }
                        string residue; //get residue
                        residue.append(1,line[17]);residue.append(1,line[18]);residue.append(1,line[19]);
                        temp_atom.residue=residue;
                        m_chains[curr_chain_ID].push_back(temp_atom);
                    }
                }
            }
        }
        file.close();
        cout<<"\nNumber of chains in protein: "<<m_chains.size()<<"\n";
        //get mesh
        m_have_mesh_atom=false;
        if(m_want_mesh_atom)
        {
            if(loadMesh(file_name.c_str())) m_have_mesh_atom=true;
        }
        m_good=true;
        calc_center_pos();//center viewing pos
    }
    //Viewing set default values
    m_rot_hori=90;m_rot_vert=90;
    m_leakage_option=0;
    m_leakage_waiting_for_response=false;
    m_marking_mode=false;
    m_black_background=true;
    m_mesh_transparency=false;
    m_anaglyph=false;
    m_interlaced=false;
    m_switch_eyes=true;
    m_eye_separation=2;
    m_convergence=30;
    m_draging=m_moving=m_drag_zooming=false;
    m_zoom=-50; //Distance between protein and viewing position
    m_window_height=m_screen_height=480;
    m_window_length=m_screen_length=640;
    m_axis_roll=0;
    init_arcball();
    m_view_mode=view_locked_roll;
    m_show_atoms=m_show_grid_protein=m_show_channel_dist=
    m_show_channel_ruler=m_show_grid_void=m_show_concave_channel=
    m_show_other_atoms=m_show_grid_void_exposure=m_show_labels=false;
    m_show_chains=m_show_mesh_atom=true;
    m_channel_edges_exp_size=0.5;
    m_channel_edges_step_size=0.5;
    m_channel_mesh_see_through=false;
    //Set protein default values
    m_grid_cut_off_left=m_grid_cut_off_right=0;
    m_pGrid=0;
    m_pGrid_color=0;
    m_channel_edge_set=false;
    m_void_vol_set=false;
    m_axis_roll_progress=0;
    m_grid_z_progress=0;
    m_void_vol_progress=0;
    m_key_delay=0;
    m_side_shift_view=0;
    m_calc_done=false;
    m_have_convex_channel_lengths=false;
    m_have_mesh_channel=false;
    m_have_mesh_normals=false;
    m_have_real_mesh_normals=false;
    m_have_mesh_channel_normals=false;
    m_have_real_mesh_channel_normals=false;
    m_mesh_normal_progress=0;
    m_mesh_normal_vertex_progress=0;
    m_solvent_radius=1.4;
    m_grid_gap_size=1.0;
    m_grid_updated=false;
    m_channel_is_leaking=false;
    m_tried_to_get_channel_mesh=false;
    m_marked_block[0]=m_marked_block[1]=m_marked_block[2]=0;
    m_channel_start_set=m_channel_end_set=m_channel_center_set=m_have_grid_void_color=false;
    m_show_grid_void_exposure_level=0;
    //Misc
}

protein::operator bool()
{
    return m_good;
}

//Calculations

void protein::calc_center_pos(void)
{
    //Calculate the center of the protein
    m_max_pos[0]=m_atoms[0].x; m_max_pos[1]=m_atoms[0].y; m_max_pos[2]=m_atoms[0].z;
    m_min_pos[0]=m_atoms[0].x; m_min_pos[1]=m_atoms[0].y; m_min_pos[2]=m_atoms[0].z;
    //get extreme values
    for (int atom=0;atom<(int)m_atoms.size();atom++)
    {
        if(m_atoms[atom].x>m_max_pos[0]) m_max_pos[0]=m_atoms[atom].x;
        if(m_atoms[atom].x<m_min_pos[0]) m_min_pos[0]=m_atoms[atom].x;
        if(m_atoms[atom].y>m_max_pos[1]) m_max_pos[1]=m_atoms[atom].y;
        if(m_atoms[atom].y<m_min_pos[1]) m_min_pos[1]=m_atoms[atom].y;
        if(m_atoms[atom].z>m_max_pos[2]) m_max_pos[2]=m_atoms[atom].z;
        if(m_atoms[atom].z<m_min_pos[2]) m_min_pos[2]=m_atoms[atom].z;
    }
    //calc center
    m_center_pos[0]=(m_max_pos[0]+m_min_pos[0])/2;
    m_center_pos[1]=(m_max_pos[1]+m_min_pos[1])/2;
    m_center_pos[2]=(m_max_pos[2]+m_min_pos[2])/2;
}

void protein::calc_grid(void)
{
    //create grid depending on protein size
    if(!m_pGrid) //initallization
    {
        m_show_grid_protein=true;
        m_protein_volume=0;
        //Set grid size (one block = 1 Å^3)
        m_grid_size[0]=int(m_max_pos[0]-m_min_pos[0]);
        m_grid_size[1]=int(m_max_pos[1]-m_min_pos[1]);
        m_grid_size[2]=int(m_max_pos[2]-m_min_pos[2]);
        //change to user specific grid box size
        m_grid_size[0]=((float)m_grid_size[0]+8)/m_grid_gap_size; // +10 to make grid larger than protein
        m_grid_size[1]=((float)m_grid_size[1]+8)/m_grid_gap_size; // +10 to make grid larger than protein
        m_grid_size[2]=((float)m_grid_size[2])/m_grid_gap_size;   // Will not be larger, causes leaks
        //allocate memory for 3D grid array
        m_pGrid=new bool**[m_grid_size[0]];
        for(int i=0;i<m_grid_size[0];i++)
        {
            m_pGrid[i] = new bool*[m_grid_size[1]];
            for(int j=0;j<m_grid_size[1];j++) m_pGrid[i][j] = new bool[m_grid_size[2]];
        }
        //output
        cout<<"\nSize of grid: X = "<<m_grid_size[0]<<", Y = "<<m_grid_size[1]<<", Z = "<<m_grid_size[2]<<endl;

        m_grid_pos=vec(m_min_pos[0]-4,
                       m_min_pos[1]-4,
                       m_min_pos[2]); //save grid position (offset 5 blocks)

        //reset grid
        for(int x=0;x<m_grid_size[0];x++)
        for(int y=0;y<m_grid_size[1];y++)
        for(int z=0;z<m_grid_size[2];z++) m_pGrid[x][y][z]=false;
    }

    //put protein on grid
    int z=m_grid_z_progress;
    for(int x=0;x<m_grid_size[0];x++)
    {
        for(int y=0;y<m_grid_size[1];y++)
        {
            if(collision_water(m_grid_pos.x+(float)x*m_grid_gap_size,
                               m_grid_pos.y+(float)y*m_grid_gap_size,
                               m_grid_pos.z+(float)z*m_grid_gap_size))//check for collision
            {
                m_pGrid[x][y][z]=true;
                m_protein_volume++;
            }
        }
    }
    if(++m_grid_z_progress>=m_grid_size[2]) //done
    {
        //cut-off fill in
        for(int i=0;i<m_grid_cut_off_left;i++)
        {
            for(int x=0;x<m_grid_size[0];x++)
            {
                for(int y=0;y<m_grid_size[1];y++)
                {
                    m_pGrid[x][y][i]=true;
                }
            }
        }
        for(int i=m_grid_size[2]-1;i>m_grid_size[2]-m_grid_cut_off_right;i--)
        {
            for(int x=0;x<m_grid_size[0];x++)
            {
                for(int y=0;y<m_grid_size[1];y++)
                {
                    m_pGrid[x][y][i]=true;
                }
            }
        }
        //fill in marked blocks
        for(int i=0;i<(int)m_void_blocks_marked.size();i++)
        {
            m_pGrid[m_void_blocks_marked[i].block_pos[0]]
                   [m_void_blocks_marked[i].block_pos[1]]
                   [m_void_blocks_marked[i].block_pos[2]]=true;
        }

        cout<<"Protein grid calculated\n";
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 11); //set text color
        float volume_factor=(1/m_grid_gap_size)*(1/m_grid_gap_size)*(1/m_grid_gap_size);
        float volume=(float)m_protein_volume/volume_factor;
        cout<<"Grid Volume of Protein is "<<volume<<" A^3"<<endl;
        SetConsoleTextAttribute(hConsole, 15); //set text color
        m_show_grid_protein=false;
        if(m_want_defined_channel)
        {
            SetConsoleTextAttribute(hConsole, 14); //set text color
            cout<<"You have selected to define the position and axis of the channel.\n"
                  "Controls: X-Axis (RED), Y-Axis (GREEN), Z-Axis (BLUE)\n"
                  "[Numpad 4] - Move marker in positive X direction\n"
                  "[Numpad 1] - Move marker in negative X direction\n"
                  "[Numpad 5] - Move marker in positive Y direction\n"
                  "[Numpad 2] - Move marker in negative Y direction\n"
                  "[Numpad 6] - Move marker in positive Z direction\n"
                  "[Numpad 3] - Move marker in negative Z direction\n"
                  "[Space]    - Mark/unmark selected grid block (red)\n"
                  "[Enter]    - When marking is complete\n"
                  "\nNow select the starting point of the channel\n";
            SetConsoleTextAttribute(hConsole, 15); //set text color
            SetForegroundWindow(m_hWnd_gui);
            m_show_grid_void=true;
            m_marked_block[0]=(float)m_grid_size[0]/2;
            m_marked_block[1]=(float)m_grid_size[1]/2;
            m_marked_block[2]=(float)m_grid_size[2]/2;
        }
        m_grid_updated=true;
    }
}

bool protein::calc_void_volumes_alt(bool exposure_test)
{

    //Starting with the grid's center block, expansions are made in all dimensions.
    //If the new block is FALSE (not collided with protein) this block will be
    //added to the block that are part of the center channel. From this new block
    //further expansions are made in all dimensions until no more free block can be found.

    if(m_void_blocks.empty())//init
    {
        SetForegroundWindow(m_hWnd_gui);
        cout<<"Measuring Grid Volume of Protein Channel... (Press C to Cancel)\n";
        m_show_grid_void=true;
        m_void_vol_progress=0;
        //start with center
        st_block_pos block_center;
        if(m_want_defined_channel&&m_channel_defined)
        {
            block_center.block_pos[0]=m_channel_center_pos.x;
            block_center.block_pos[1]=m_channel_center_pos.y;
            block_center.block_pos[2]=m_channel_center_pos.z;
        }
        else
        {
            block_center.block_pos[0]=m_grid_size[0]/2;
            block_center.block_pos[1]=m_grid_size[1]/2;
            block_center.block_pos[2]=m_grid_size[2]/2;
        }
        m_void_blocks_exp.push_back(block_center);
        //begin expansion (center block must be free!)
        bool found_cleft=true;
        if(m_pGrid[block_center.block_pos[0]][block_center.block_pos[1]][block_center.block_pos[2]])
        {
            found_cleft=false;
            //center is not void
            cout<<"Center of the channel is blocked, will try to find nearby channel\n";
            for(int exp=1;exp<1000;exp++)//will expand sideways to find channel
            {
                //x+
                if(block_center.block_pos[0]+exp<m_grid_size[0])
                {
                    if(!m_pGrid[block_center.block_pos[0]+exp][block_center.block_pos[1]][block_center.block_pos[2]])
                    {
                        found_cleft=true;
                        block_center.block_pos[0]+=exp;
                        break;
                    }
                }
                //x-
                if(block_center.block_pos[0]-exp>=0)
                {
                    if(!m_pGrid[block_center.block_pos[0]-exp][block_center.block_pos[1]][block_center.block_pos[2]])
                    {
                        found_cleft=true;
                        block_center.block_pos[0]-=exp;
                        break;
                    }
                }
                //y+
                if(block_center.block_pos[1]+exp<m_grid_size[1])
                {
                    if(!m_pGrid[block_center.block_pos[0]][block_center.block_pos[1]+exp][block_center.block_pos[2]])
                    {
                        found_cleft=true;
                        block_center.block_pos[1]+=exp;
                        break;
                    }
                }
                //y-
                if(block_center.block_pos[1]-exp>=0)
                {
                    if(!m_pGrid[block_center.block_pos[0]][block_center.block_pos[1]-exp][block_center.block_pos[2]])
                    {
                        found_cleft=true;
                        block_center.block_pos[1]-=exp;
                        break;
                    }
                }
                //z+
                if(block_center.block_pos[2]+exp<m_grid_size[2])
                {
                    if(!m_pGrid[block_center.block_pos[0]][block_center.block_pos[1]][block_center.block_pos[2]+exp])
                    {
                        found_cleft=true;
                        block_center.block_pos[2]+=exp;
                        break;
                    }
                }
                //z-
                if(block_center.block_pos[2]-exp>=0)
                {
                    if(!m_pGrid[block_center.block_pos[0]][block_center.block_pos[1]][block_center.block_pos[2]-exp])
                    {
                        found_cleft=true;
                        block_center.block_pos[2]-=exp;
                        break;
                    }
                }
            }
        }
        if(found_cleft) m_void_blocks.push_back(block_center);
        else//give up
        {
            HANDLE hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
            SetConsoleTextAttribute(hConsole, 12); //set RED text color
            cout<<"Could not find a cleft near center of the channel.\nPlease redefine location of the channel\n";
            SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
            m_want_void_vol=false;
            m_want_mesh_channel=false;
            return false;
        }
    }

    int index=m_void_vol_progress; //index of curret position in void_blocks_exp
    //origin
    int x_pos,y_pos,z_pos;
    x_pos=m_void_blocks_exp[index].block_pos[0];
    y_pos=m_void_blocks_exp[index].block_pos[1];
    z_pos=m_void_blocks_exp[index].block_pos[2];
    //leak test (that grid walls is reached indicates that the channel is leaking)
    if(!m_channel_is_leaking)
    {
        if(x_pos==0 || x_pos==m_grid_size[0] ||
           y_pos==0 || y_pos==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
        }
    }
    if(m_channel_is_leaking)
    {
        if(!m_leakage_waiting_for_response&&m_leakage_option==0)//show once
        {
            //SetForegroundWindow(m_hWnd_console);
            HANDLE hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
            SetConsoleTextAttribute(hConsole, 12); //set text color
            cout<<"\n***************************\nSolvent can pass through the Channel's sides!\n***************************\n";
            SetConsoleTextAttribute(hConsole, 14); //set text color
            cout<<
            "To continue the leakage must be stopped by one of the following options:\n"
            "[1] - Change the radius of the solvet compound\n"
            "[2] - Define position and direction of the Channel/Cleft\n"
            "[3] - Manually select grid blocks to block the pathway of the leakage\n\n"
            "Choose one option by pressing the corresponding number on your keyboard\nwhile the protein display window is marked.\n";
            m_leakage_waiting_for_response=true;
        }
        if(!m_leakage_waiting_for_response&&m_leakage_option!=0)//enter when option is choosen
        {
            switch(m_leakage_option)
            {
                case 1://Change solvent size
                {
                    HANDLE hConsole;
                    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                    SetForegroundWindow(m_hWnd_console);//focus console
                    SetConsoleTextAttribute(hConsole, 14); //set text color
                    cout<<"\nPlease enter radius of solvent (default is 1.4 Angstrom): ";
                    SetConsoleTextAttribute(hConsole, 15); //set text color
                    //input
                    string sInput;
                    getline(cin,sInput);//get new radius
                    cout<<endl;
                    float solvent_radius=atof(sInput.c_str());
                    if(solvent_radius>0) m_solvent_radius=solvent_radius;//must be positive
                    else return false;
                    //clean-up
                    m_void_blocks.clear();
                    m_void_blocks_exp.clear();
                    m_void_blocks_new.clear();

                    //reset grid
                    for(int x=0;x<m_grid_size[0];x++)
                    for(int y=0;y<m_grid_size[1];y++)
                    for(int z=0;z<m_grid_size[2];z++) m_pGrid[x][y][z]=false;

                    m_grid_updated=false;
                    m_grid_z_progress=0;
                    m_void_vol_progress=0;
                    m_protein_volume=0;

                    m_show_grid_protein=true;
                    m_channel_is_leaking=false;
                    m_leakage_option=0;
                    SetForegroundWindow(m_hWnd_gui);
                    return false;
                }break;
                case 2://Select channel
                {
                    m_want_defined_channel=true;
                    m_channel_defined=false;
                    HANDLE hConsole;
                    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                    //SetForegroundWindow(m_hWnd_console);//focus console
                    SetConsoleTextAttribute(hConsole, 14); //set text color
                    cout<<"You have selected to define the position and direction of the channel.\n"
                          "Controls: X-Axis (RED), Y-Axis (GREEN), Z-Axis (BLUE)\n"
                          "[Numpad 4] - Move marker in positive X direction\n"
                          "[Numpad 1] - Move marker in negative X direction\n"
                          "[Numpad 5] - Move marker in positive Y direction\n"
                          "[Numpad 2] - Move marker in negative Y direction\n"
                          "[Numpad 6] - Move marker in positive Z direction\n"
                          "[Numpad 3] - Move marker in negative Z direction\n"
                          "[Space]    - Mark selected grid block\n"
                          "[Del]      - Remove selection\n"
                          "[Enter]    - When marking is complete\n"
                          "\nNow select the starting point of the channel\n";
                    SetConsoleTextAttribute(hConsole, 15); //set text color
                    m_show_grid_void=true;
                    m_marked_block[0]=(float)m_grid_size[0]/2;
                    m_marked_block[1]=(float)m_grid_size[1]/2;
                    m_marked_block[2]=(float)m_grid_size[2]/2;

                    //clean-up
                    m_void_blocks.clear();
                    m_void_blocks_exp.clear();
                    m_void_blocks_new.clear();
                    m_void_vol_progress=0;

                    m_channel_is_leaking=false;
                    m_leakage_option=0;
                    return false;
                }break;
                case 3://Block manually
                {
                    if(!m_marking_mode)
                    {
                        m_marking_mode=true;
                        HANDLE hConsole;
                        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                        SetConsoleTextAttribute(hConsole, 14); //set text color
                        cout<<"\nYou can now manually choose grid blocks that not will be part of the channel\n"
                              "Controls:\n"
                              "[Numpad 4] - Move marker in positive X direction\n"
                              "[Numpad 1] - Move marker in negative X direction\n"
                              "[Numpad 5] - Move marker in positive Y direction\n"
                              "[Numpad 2] - Move marker in negative Y direction\n"
                              "[Numpad 6] - Move marker in positive Z direction\n"
                              "[Numpad 3] - Move marker in negative Z direction\n"
                              "[Space]    - Mark/unmark selected grid block (red)\n"
                              "[Enter]    - When marking is complete\n";
                        SetConsoleTextAttribute(hConsole, 15); //set text color
                        //set marked block to center
                        m_marked_block[0]=m_grid_size[0]/2;
                        m_marked_block[1]=m_grid_size[1]/2;
                        m_marked_block[2]=m_grid_size[2]/2;
                    }
                }break;
                case 4://Manual blocking done
                {
                    m_marking_mode=false;
                    //clean-up
                    m_void_blocks.clear();
                    m_void_blocks_exp.clear();
                    m_void_blocks_new.clear();
                    //reset grid
                    for(int x=0;x<m_grid_size[0];x++)
                    for(int y=0;y<m_grid_size[1];y++)
                    for(int z=0;z<m_grid_size[2];z++) m_pGrid[x][y][z]=false;
                    m_grid_updated=false;
                    m_grid_z_progress=0;
                    m_void_vol_progress=0;
                    m_protein_volume=0;

                    m_show_grid_protein=true;
                    m_channel_is_leaking=false;
                    m_leakage_option=0;
                    SetForegroundWindow(m_hWnd_gui);
                    return false;
                }
            }
        }
        return false;
    }
    //expansions
    st_block_pos block_temp;
    int cut_box=0;//0 = OFF
    //x+
    for(int exp=1; x_pos+exp<m_grid_size[0]-cut_box; exp++)
    {
        if(m_pGrid[x_pos+exp][y_pos][z_pos]) break; //hit the wall
        //leak test
        if(x_pos+exp==0 || x_pos+exp==m_grid_size[0] || y_pos==0 || y_pos==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos+exp][y_pos][z_pos]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos+exp;
        block_temp.block_pos[1]=y_pos;
        block_temp.block_pos[2]=z_pos;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);//block is part of the channel
    }
    //x-
    for(int exp=1; x_pos-exp>=cut_box; exp++)
    {
        if(m_pGrid[x_pos-exp][y_pos][z_pos]) break; //hit the wall
        //leak test
        if(x_pos-exp==0 || x_pos-exp==m_grid_size[0] || y_pos==0 || y_pos==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos-exp][y_pos][z_pos]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos-exp;
        block_temp.block_pos[1]=y_pos;
        block_temp.block_pos[2]=z_pos;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);
    }
    //y+
    for(int exp=1; y_pos+exp<m_grid_size[1]-cut_box; exp++)
    {
        if(m_pGrid[x_pos][y_pos+exp][z_pos]) break;
        //leak test
        if(x_pos==0 || x_pos==m_grid_size[0] || y_pos+exp==0 || y_pos+exp==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos][y_pos+exp][z_pos]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos;
        block_temp.block_pos[1]=y_pos+exp;
        block_temp.block_pos[2]=z_pos;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);
    }
    //y-
    for(int exp=1; y_pos-exp>=cut_box; exp++)
    {
        if(m_pGrid[x_pos][y_pos-exp][z_pos]) break;
        //leak test
        if(x_pos==0 || x_pos==m_grid_size[0] || y_pos-exp==0 || y_pos-exp==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos][y_pos-exp][z_pos]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos;
        block_temp.block_pos[1]=y_pos-exp;
        block_temp.block_pos[2]=z_pos;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);
    }
    //z+
    for(int exp=1; z_pos+exp<m_grid_size[2]-cut_box; exp++)
    {
        if(m_pGrid[x_pos][y_pos][z_pos+exp]) break;
        //leak test
        if(x_pos==0 || x_pos==m_grid_size[0] || y_pos==0 || y_pos==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos][y_pos][z_pos+exp]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos;
        block_temp.block_pos[1]=y_pos;
        block_temp.block_pos[2]=z_pos+exp;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);
    }
    //z-
    for(int exp=1; z_pos-exp>=cut_box; exp++)
    {
        if(m_pGrid[x_pos][y_pos][z_pos-exp]) break;
        //leak test
        if(x_pos==0 || x_pos==m_grid_size[0] || y_pos==0 || y_pos==m_grid_size[1])
        {
            m_channel_is_leaking=true;
            m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
            m_void_blocks_new.clear();
            return false;
        }
        //exposure level test
        if(exposure_test&&m_pGrid_color[x_pos][y_pos][z_pos-exp]>m_show_grid_void_exposure_level) break;

        block_temp.block_pos[0]=x_pos;
        block_temp.block_pos[1]=y_pos;
        block_temp.block_pos[2]=z_pos-exp;
        //if channel is user defined test if block is inside planes
        if(m_want_defined_channel&&m_channel_defined) if(!is_block_inside_planes(block_temp)) break;
        bool already_stored=false;
        for(int block_nr=0;block_nr<(int)m_void_blocks.size();block_nr++)
        {
            if(block_temp.block_pos[0]==m_void_blocks[block_nr].block_pos[0] &&
               block_temp.block_pos[1]==m_void_blocks[block_nr].block_pos[1] &&
               block_temp.block_pos[2]==m_void_blocks[block_nr].block_pos[2]) //block already in store
            {
                already_stored=true;
                break;
            }
        }
        if(!already_stored) //also check that block is now already in new blocks vector
        {
            for(int block_nr=0;block_nr<(int)m_void_blocks_new.size();block_nr++)
            {
                if(block_temp.block_pos[0]==m_void_blocks_new[block_nr].block_pos[0] &&
                   block_temp.block_pos[1]==m_void_blocks_new[block_nr].block_pos[1] &&
                   block_temp.block_pos[2]==m_void_blocks_new[block_nr].block_pos[2]) //block already in store
                {
                    already_stored=true;
                    break;
                }
            }
        }
        if(already_stored)continue;
        else m_void_blocks_new.push_back(block_temp);
    }
    //if new vector is empty=done
    if(m_void_blocks_new.empty())
    {
        if(m_void_vol_progress==(int)m_void_blocks_exp.size()-1)//only don if all blocks in exp vector has been checked (on the last element)
        {
            //output
            HANDLE hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
            SetConsoleTextAttribute(hConsole, 11); //set text color
            float volume_factor=(1/m_grid_gap_size)*(1/m_grid_gap_size)*(1/m_grid_gap_size);
            float volume=(float)m_void_blocks.size()/volume_factor;
            cout<<"Grid Volume of Protein Channel is "<<volume<<" A^3\n";
            SetConsoleTextAttribute(hConsole, 15); //set text color
            //clean-up
            m_void_blocks_exp.clear();
            m_void_blocks_new.clear();
            m_void_blocks_marked.clear();
            m_void_vol_set=true;
            m_show_grid_void=false;
        }
    }
    else if(m_void_vol_progress==(int)m_void_blocks_exp.size()-1) //on the last element in m_void_blocks_exp (but not done), put in new blocks from new vector
    {
        //empty exp vector
        m_void_blocks_exp.clear();
        //move new vector to void and exp vector
        m_void_blocks_exp.insert(m_void_blocks_exp.begin(),m_void_blocks_new.begin(),m_void_blocks_new.end());
        m_void_blocks.insert(m_void_blocks.end(),m_void_blocks_new.begin(),m_void_blocks_new.end());
        //empty new vector
        m_void_blocks_new.clear();
        //reset the indexing counter
        m_void_vol_progress=-1;
    }

    m_void_vol_progress++; //increment index counter
    return true;
}

bool protein::calc_grid_void_color(void)
{
    //Test how exposed grid blocks are
    if(!m_pGrid_color) //initallization
    {
        //allocate memory for 3D grid array
        m_pGrid_color=new int**[m_grid_size[0]];
        for(int i=0;i<m_grid_size[0];i++)
        {
            m_pGrid_color[i] = new int*[m_grid_size[1]];
            for(int j=0;j<m_grid_size[1];j++) m_pGrid_color[i][j] = new int[m_grid_size[2]];
        }
        //reset grid color
        for(int x=0;x<m_grid_size[0];x++)
        for(int y=0;y<m_grid_size[1];y++)
        for(int z=0;z<m_grid_size[2];z++) m_pGrid_color[x][y][z]=6;

        m_grid_z_progress=1;//0 is part of wall
    }
    int z=m_grid_z_progress;
    for(int x=1;x<m_grid_size[0]-1;x++)
    {
        for(int y=1;y<m_grid_size[1]-1;y++)
        {
            if(m_pGrid[x][y][z])//part of protein
            {
                m_pGrid_color[x][y][z]=7;
                continue;
            }
            //expansion test
            int num_of_walls_hit=0;
            //x+
            for(int exp=1;x+exp<m_grid_size[0];exp++)
            {
                if(m_pGrid[x+exp][y][z]) break;//protein collision
                if(x+exp==m_grid_size[0]-1) {num_of_walls_hit++;break;}//wall collision
            }
            //x-
            for(int exp=1;x-exp>=0;exp++)
            {
                if(m_pGrid[x-exp][y][z]) break;//protein collision
                if(x-exp==0)                {num_of_walls_hit++;break;}//wall collision
            }
            //y+
            for(int exp=1;y+exp<m_grid_size[1];exp++)
            {
                if(m_pGrid[x][y+exp][z]) break;//protein collision
                if(y+exp==m_grid_size[1]-1) {num_of_walls_hit++;break;}//wall collision
            }
            //y-
            for(int exp=1;y-exp>=0;exp++)
            {
                if(m_pGrid[x][y-exp][z]) break;//protein collision
                if(y-exp==0)                {num_of_walls_hit++;break;}//wall collision
            }
            //z+
            for(int exp=1;z+exp<m_grid_size[2];exp++)
            {
                if(m_pGrid[x][y][z+exp]) break;//protein collision
                if(z+exp==m_grid_size[2]-1) {num_of_walls_hit++;break;}//wall collision
            }
            //z-
            for(int exp=1;z-exp>=0;exp++)
            {
                if(m_pGrid[x][y][z-exp]) break;//protein collision
                if(z-exp==0)                {num_of_walls_hit++;break;}//wall collision
            }
            //assing value for exposure
            m_pGrid_color[x][y][z]=num_of_walls_hit;
        }
    }

    if(++m_grid_z_progress>=m_grid_size[2]) //done
    {
        cout<<"Level of Exposure for Clefts calculated\n";
        m_have_grid_void_color=true;
        m_show_grid_void_exposure=false;
    }
    return true;
}

bool protein::calc_channel_edges(void)
{
    //From the axis along the channel expansions are made to the sides (perpendicular to the axis)
    //until the edge of the channel is struck. The final point is stored to be used in the channel mesh.
    if(!m_show_channel_ruler) m_show_channel_ruler=true;//to view measurment progress

    float step_size_axis=m_channel_edges_step_size;//in Å
    float step_size_sides=m_channel_edges_exp_size;//in Å
    float max_exp=20;//to avoid neverending expansions of a hole is found
    int roll=m_axis_roll=m_axis_roll_progress;//get current angle

    if(m_want_defined_channel&&m_channel_defined)
    {
        //Init
        if(m_channel_edges_up.empty())
        {
            m_axis_roll_progress=0;
            //get number of steps in measurments
            float dx=(m_channel_end_pos.x-m_channel_start_pos.x)*m_grid_gap_size;
            float dy=(m_channel_end_pos.y-m_channel_start_pos.y)*m_grid_gap_size;
            float dz=(m_channel_end_pos.z-m_channel_start_pos.z)*m_grid_gap_size;
            float channel_length=sqrt(dx*dx+dy*dy+dz*dz);
            m_channel_sidesteps=channel_length/m_channel_edges_step_size;
            m_channel_expsteps=channel_length/m_channel_edges_exp_size;
        }

        //create channel vector
        vec axis_norm=m_channel_start_pos-m_channel_end_pos;
        axis_norm=axis_norm.unit();
        //create vector for cross product
        vec vec_rand(-1,0,0);
        if(axis_norm.equals(vec_rand)||axis_norm.equals(vec(1,0,0))) //will give 2 parallell vectors
        {
            vec_rand=vec(0,0,1);//will not be parallell
        }
        //Move point around axis for the current progress angle (roll)
        vec rotated_vec=rotatePointAboutLine(vec_rand,float(roll)*_piover180,vec(0,0,0),axis_norm);
        //create perpendicular vector
        vec cross_vec=m_channel_vector^rotated_vec;
        cross_vec=cross_vec.unit();

        //direction done
        float xr=cross_vec.x*step_size_sides;
        float yr=cross_vec.y*step_size_sides;
        float zr=cross_vec.z*step_size_sides;

        //up measurment
        m_channel_edges_up.push_back(vector<vec>());
        for(int axis_step=0;axis_step<m_channel_sidesteps;axis_step++)
        {
            //start pos
            float xs=float(m_channel_start_pos.x)*m_grid_gap_size+m_grid_pos.x+m_channel_vector.x*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size; //due to vector is half grid size
            float ys=float(m_channel_start_pos.y)*m_grid_gap_size+m_grid_pos.y+m_channel_vector.y*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size;
            float zs=float(m_channel_start_pos.z)*m_grid_gap_size+m_grid_pos.z+m_channel_vector.z*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size;

            for(float exp_step=0;exp_step<=m_channel_expsteps;exp_step+=1)
            {
                if(collision_water(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr))
                {
                    m_channel_edges_up[roll].push_back(vec(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr));
                    break;
                }
                if(exp_step>m_channel_expsteps-1) //maximum expansion reached
                {
                    m_channel_edges_up[roll].push_back(vec(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr));
                    break;
                }
            }
        }
        //down measurment (+180°)
        xr=-xr;
        yr=-yr;
        zr=-zr;
        m_channel_edges_down.push_back(vector<vec>());
        for(int axis_step=0;axis_step<m_channel_sidesteps;axis_step++)
        {
            //start pos
            float xs=float(m_channel_start_pos.x)*m_grid_gap_size+m_grid_pos.x+m_channel_vector.x*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size; //due to vector is half grid size
            float ys=float(m_channel_start_pos.y)*m_grid_gap_size+m_grid_pos.y+m_channel_vector.y*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size;
            float zs=float(m_channel_start_pos.z)*m_grid_gap_size+m_grid_pos.z+m_channel_vector.z*(float)axis_step*2/m_grid_gap_size*m_channel_edges_step_size;

            for(float exp_step=0;exp_step<=m_channel_expsteps;exp_step+=1)
            {
                if(collision_water(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr))
                {
                    m_channel_edges_down[roll].push_back(vec(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr));
                    break;
                }
                if(exp_step>m_channel_expsteps-1) //maximum expansion reached
                {
                    m_channel_edges_down[roll].push_back(vec(xs+xr*exp_step,ys+yr*exp_step,zs+exp_step*zr));
                    break;
                }
            }
        }
    }
    else//channel located along Z-axis
    {
        //Init
        if(m_channel_edges_up.empty())
        {
            m_axis_roll_progress=0;
            //get number of steps in measurments
            float dx=(m_channel_end_pos.x-m_channel_start_pos.x)*m_grid_gap_size;
            float dy=(m_channel_end_pos.y-m_channel_start_pos.y)*m_grid_gap_size;
            float dz=(m_channel_end_pos.z-m_channel_start_pos.z)*m_grid_gap_size;
            float channel_length=sqrt(dx*dx+dy*dy+dz*dz);
            m_channel_sidesteps=channel_length/m_channel_edges_step_size;
            m_channel_expsteps=channel_length/m_channel_edges_exp_size;
            //if(!m_channel_defined)
        }

        float xr=sinf(float(roll)*_piover180);
        float yr=cosf(float(roll)*_piover180);
        //up measurment
        m_channel_edges_up.push_back(vector<vec>());
        for(float axis_pos=m_min_pos[2];axis_pos<m_max_pos[2];axis_pos+=step_size_axis)
        {
            for(float exp=0;exp<max_exp;exp+=step_size_sides)
            {

                if(collision_water(xr*exp,yr*exp,axis_pos))
                {
                    m_channel_edges_up[roll].push_back(vec(xr*exp,yr*exp,axis_pos));
                    break;
                }
                if(exp>max_exp-step_size_sides*2) //maximum expansion reached
                {
                    m_channel_edges_up[roll].push_back(vec(xr*exp,yr*exp,axis_pos));
                    break;
                }
            }
        }
        xr=-xr;
        yr=-yr;
        //down measurment (+180°)
        m_channel_edges_down.push_back(vector<vec>());
        for(float axis_pos=m_min_pos[2];axis_pos<m_max_pos[2];axis_pos+=step_size_axis)
        {
            for(float exp=0;exp<max_exp;exp+=step_size_sides)
            {

                if(collision_water(xr*exp,yr*exp,axis_pos))
                {
                    m_channel_edges_down[roll].push_back(vec(xr*exp,yr*exp,axis_pos));
                    break;
                }
                if(exp>max_exp-step_size_sides*2) //maximum expansion reached
                {
                    m_channel_edges_down[roll].push_back(vec(xr*exp,yr*exp,axis_pos));
                    break;
                }
            }
        }
    }
    //check if done
    if(++m_axis_roll_progress>=180)
    {
        m_channel_edge_set=true;
        m_show_channel_ruler=false;
        m_show_channel_dist=true;
        //merge up and down vector
        for(int segment=0;segment<(int)m_channel_edges_down.size();segment++)
        {
            m_channel_edges_up.push_back(m_channel_edges_down[segment]);
        }
        //push backs for normals
        for(int iRoll=0;iRoll<(int)m_channel_edges_up.size();iRoll++)
        {
            m_channel_edges_normals.push_back(vector<vec>());//push back normal vector for comming normals
            for(int segment=0;segment<(int)m_channel_edges_normals[iRoll].size();segment++)
            {
                m_channel_edges_normals[iRoll].push_back(vec(0,0,0));//add empty normal
            }
        }
        m_channel_edges_down.clear();
        cout<<"Channel Diameters Measured\n";
    }
    return true;
}

bool protein::calc_mesh_normals(void)
{
    //calc common normals for all triangles
    int max_rounds=10000;
    for(int round=0;round<max_rounds;round++)
    {
        int triangle=m_mesh_normal_progress;
        //store vertex to triangle info in map
        m_vertex_triangle_map_atom[m_mesh_atom_triangles[triangle].first].push_back(triangle);
        m_vertex_triangle_map_atom[m_mesh_atom_triangles[triangle].second].push_back(triangle);
        m_vertex_triangle_map_atom[m_mesh_atom_triangles[triangle].third].push_back(triangle);
        //calculate normal for first and assign to second and third
        vec center(m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].first].x,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].first].y,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].first].z);
        vec vec1  (m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].second].x,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].second].y,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].second].z);
        vec vec2  (m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].third].x,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].third].y,
                   m_mesh_atom_vertices[m_mesh_atom_triangles[triangle].third].z);
        vec1=vec1-center;
        vec2=vec2-center;
        vec normal=vec1^vec2;
        normal=normal.unit();//normalise
        //assigning
        m_mesh_atom_triangles[triangle].first_normal=normal;
        m_mesh_atom_triangles[triangle].second_normal=normal;
        m_mesh_atom_triangles[triangle].third_normal=normal;
        //check if done
        if(++m_mesh_normal_progress>=(int)m_mesh_atom_triangles.size())
        {
            m_have_mesh_normals=true;
            cout<<"Common mesh normals calculated\n";
            return true;
        }
    }
    return true;
}

bool protein::calc_mesh_channel_normals(void)
{
    //calc common normals for all triangles
    int max_rounds=1000;
    for(int round=0;round<max_rounds;round++)
    {
        int triangle=m_mesh_normal_progress;
        //store vertex to triangle info in map
        m_vertex_triangle_map_channel[m_mesh_channel_triangles[triangle].first].push_back(triangle);
        m_vertex_triangle_map_channel[m_mesh_channel_triangles[triangle].second].push_back(triangle);
        m_vertex_triangle_map_channel[m_mesh_channel_triangles[triangle].third].push_back(triangle);
        //calculate normal for first and assign to second and third
        vec center(m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].first].x,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].first].y,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].first].z);
        vec vec1  (m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].second].x,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].second].y,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].second].z);
        vec vec2  (m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].third].x,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].third].y,
                   m_mesh_channel_vertices[m_mesh_channel_triangles[triangle].third].z);
        vec1=vec1-center;
        vec2=vec2-center;
        vec normal=vec1^vec2;
        normal=normal.unit();//normalise
        //assigning
        m_mesh_channel_triangles[triangle].first_normal=normal;
        m_mesh_channel_triangles[triangle].second_normal=normal;
        m_mesh_channel_triangles[triangle].third_normal=normal;
        //check if done
        if(++m_mesh_normal_progress>=(int)m_mesh_channel_triangles.size())
        {
            m_have_mesh_channel_normals=true;
            cout<<"Common channel mesh normals calculated\n";
            return true;
        }
    }
    return true;
}

bool protein::calc_real_mesh_normals(void)
{
    //calc specific normals for all vertices
    int max_rounds=10000;//number of normals calculated per cycle
    for(int round=0;round<max_rounds;round++)
    {
        int vertex=m_mesh_normal_vertex_progress;
        vector<vec> normals;
        //get normals for current vertex
        for(int index=0;index<(int)m_vertex_triangle_map_atom[vertex].size();index++)//go through all triangles involved
        {
            int triangle=m_vertex_triangle_map_atom[vertex][index];
            //find if vertex belongs to first,second or third vertex in triangle
            if(m_mesh_atom_triangles[triangle].first==vertex)
            {
                normals.push_back(m_mesh_atom_triangles[triangle].first_normal);
            }
            if(m_mesh_atom_triangles[triangle].second==vertex)
            {
                normals.push_back(m_mesh_atom_triangles[triangle].second_normal);
            }
            if(m_mesh_atom_triangles[triangle].third==vertex)
            {
                normals.push_back(m_mesh_atom_triangles[triangle].third_normal);
            }
        }
        vec sum(0,0,0);
        for(int normal_index=0;normal_index<(int)normals.size();normal_index++)//sum normals to average
        {
            sum=sum+normals[normal_index];
        }
        vec real_normal=sum.unit();//normalise
        //assign real normal
        for(int index=0;index<(int)m_vertex_triangle_map_atom[vertex].size();index++)//go through all triangles involved
        {
            int triangle=m_vertex_triangle_map_atom[vertex][index];
            //find if vertex belongs to first,second or third vertex in triangle
            if(m_mesh_atom_triangles[triangle].first==vertex)
            {
                m_mesh_atom_triangles[triangle].first_avg_normal=real_normal;
            }
            if(m_mesh_atom_triangles[triangle].second==vertex)
            {
                m_mesh_atom_triangles[triangle].second_avg_normal=real_normal;
            }
            if(m_mesh_atom_triangles[triangle].third==vertex)
            {
                m_mesh_atom_triangles[triangle].third_avg_normal=real_normal;
            }
        }
        if(++m_mesh_normal_vertex_progress>=(int)m_mesh_atom_vertices.size())
        {
            m_have_real_mesh_normals=true;
            m_show_mesh_atom=false;
            m_mesh_normal_vertex_progress=0;
            m_mesh_normal_progress=0;
            cout<<"Average mesh normals calculated\n";
            return true;
        }
    }
    return true;
}

bool protein::calc_real_mesh_channel_normals(void)
{
    //calc specific normals for all vertices
    int max_rounds=1000;//number of normals calculated per cycle
    for(int round=0;round<max_rounds;round++)
    {
        int vertex=m_mesh_normal_vertex_progress;
        vector<vec> normals;
        //get normals for current vertex
        for(int index=0;index<(int)m_vertex_triangle_map_channel[vertex].size();index++)//go through all triangles involved
        {
            int triangle=m_vertex_triangle_map_channel[vertex][index];
            //find if vertex belongs to first,second or third vertex in triangle
            if(m_mesh_channel_triangles[triangle].first==vertex)
            {
                normals.push_back(m_mesh_channel_triangles[triangle].first_normal);
            }
            if(m_mesh_channel_triangles[triangle].second==vertex)
            {
                normals.push_back(m_mesh_channel_triangles[triangle].second_normal);
            }
            if(m_mesh_channel_triangles[triangle].third==vertex)
            {
                normals.push_back(m_mesh_channel_triangles[triangle].third_normal);
            }
        }
        vec sum(0,0,0);
        for(int normal_index=0;normal_index<(int)normals.size();normal_index++)//sum normals to average
        {
            sum=sum+normals[normal_index];
        }
        vec real_normal=sum.unit();//normalise
        //assign real normal
        for(int index=0;index<(int)m_vertex_triangle_map_channel[vertex].size();index++)//go through all triangles involved
        {
            int triangle=m_vertex_triangle_map_channel[vertex][index];
            //find if vertex belongs to first,second or third vertex in triangle
            if(m_mesh_channel_triangles[triangle].first==vertex)
            {
                m_mesh_channel_triangles[triangle].first_avg_normal=real_normal;
            }
            if(m_mesh_channel_triangles[triangle].second==vertex)
            {
                m_mesh_channel_triangles[triangle].second_avg_normal=real_normal;
            }
            if(m_mesh_channel_triangles[triangle].third==vertex)
            {
                m_mesh_channel_triangles[triangle].third_avg_normal=real_normal;
            }
        }
        if(++m_mesh_normal_vertex_progress>=(int)m_mesh_channel_vertices.size())
        {
            m_have_real_mesh_channel_normals=true;
            m_mesh_normal_vertex_progress=0;
            cout<<"Average mesh normals calculated\n";
            return true;
        }
    }
    return true;
}

bool protein::calc_convex_channel_lengths(void)
{
    if(!m_channel_defined)
    {
        m_channel_start_pos=vec(m_grid_size[0]/2,m_grid_size[1]/2,0);
        m_channel_vector=vec(0,0,1)*(m_grid_gap_size/2);
    }
    //Measure length of the expansions inside the channel and assign a color depending on length
    float longest=0;
    for(int roll=0;roll<(int)m_channel_edges_up.size();roll++)
    {
        m_channel_edges_length.push_back(vector<float>());
        //measure lengths
        for(int segment=0;segment<(int)m_channel_edges_up[roll].size();segment++)
        {
            //measure length
            float cons=(float)segment*m_channel_edges_step_size*2/m_grid_gap_size;
            vec center=(m_channel_vector*cons)+(m_channel_start_pos*m_grid_gap_size)+m_grid_pos; //start of expansion pos
            vec edge(m_channel_edges_up[roll][segment].x,
                     m_channel_edges_up[roll][segment].y,
                     m_channel_edges_up[roll][segment].z);
            vec diff=edge-center;
            float length=diff.length();
            m_channel_edges_length[roll].push_back(length);
            if(length>longest) longest=length; //store longest value
        }
    }
    for(int roll=0;roll<(int)m_channel_edges_up.size();roll++)
    {
        m_channel_edges_color.push_back(vector<vec>());
        //assign colors
        for(int segment=0;segment<(int)m_channel_edges_up[roll].size();segment++)
        {
            float length=m_channel_edges_length[roll][segment];
            vec color(length/longest,1-length/longest,0.5);
            m_channel_edges_color[roll].push_back(color);
        }
    }
    m_have_convex_channel_lengths=true;
    cout<<"Distances across channel calculated\n";

    return true;
}

bool protein::calc_concave_channel_mesh(void)
{
    m_tried_to_get_channel_mesh=true;
    //Create a fake pdb file with water molecules in position of the channel void blocks
    ofstream file("data\\solvent_part.pdb");
    if(file==0) return false;
    file<<"RADIUS OF SOLVENT "<<m_solvent_radius<<endl;
    string start_part("SOLVENT                        ");
    for(int block=0;block<(int)m_void_blocks.size();block++)
    {
        string temp_line;
        temp_line.append(start_part);
        //x
        char cXpos[8]={'0','0','0','.','0','0','0','\0'};
        float xpos=float(m_void_blocks[block].block_pos[0])*m_grid_gap_size+m_grid_pos.x;
        if(!float_to_char(xpos,cXpos)) return false;
        temp_line.append(cXpos);
        temp_line.append(" ");//for spacing
        //y
        char cYpos[8]={'0','0','0','.','0','0','0','\0'};
        float ypos=float(m_void_blocks[block].block_pos[1])*m_grid_gap_size+m_grid_pos.y;
        if(!float_to_char(ypos,cYpos)) return false;
        temp_line.append(cYpos);
        temp_line.append(" ");//for spacing
        //z
        char cZpos[8]={'0','0','0','.','0','0','0','\0'};
        float zpos=float(m_void_blocks[block].block_pos[2])*m_grid_gap_size+m_grid_pos.z;
        if(!float_to_char(zpos,cZpos)) return false;
        temp_line.append(cZpos);

        file<<temp_line<<endl;
    }
    file.close();

    //Send the fake pdb file to MeshGenerator to get a mesh
    if(!load_mesh_channel("solvent_part.pdb")) return false;

    m_show_concave_channel=true;
    m_show_channel_dist=false;

    return true;
}

bool protein::calc_channel_mesh_vol(void)
{
    float sum=0;
    for(int i=0;i<(int)m_mesh_channel_triangles.size();i++)
    {
        sum+=volume_of_triangle(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first],
                                m_mesh_channel_vertices[m_mesh_channel_triangles[i].second],
                                m_mesh_channel_vertices[m_mesh_channel_triangles[i].third]);
    }
    m_channel_mesh_vol=sum;
    return true;
}


//Update

void protein::update_protein(void)
{
    //update calculations
    if(!m_calc_done)
    {
        if(m_have_mesh_atom && !m_have_mesh_normals) calc_mesh_normals();
        else//continue with real normals
        {
            if(m_have_mesh_atom && !m_have_real_mesh_normals) calc_real_mesh_normals();
            else//continue with grid
            {
                if(!m_grid_updated) calc_grid();
                else//continue with channel edges
                {
                    if(!m_channel_defined && m_want_defined_channel) m_marking_mode=true;
                    else
                    {
                        if(!m_have_grid_void_color) calc_grid_void_color();
                        else
                        {
                            if(!m_void_vol_set && m_want_void_vol) calc_void_volumes_alt(m_want_defined_channel);
                            else//continue measuring distances inside channel
                            {
                                if(!m_channel_edge_set && m_want_dist_channel) calc_channel_edges();
                                else//continue with normals for channel
                                {
                                    if(!m_have_convex_channel_lengths && m_want_dist_channel) calc_convex_channel_lengths();
                                    else//continue with concave channel mesh
                                    {
                                        if(!m_tried_to_get_channel_mesh && m_want_mesh_channel) //one try
                                        {
                                            if(calc_concave_channel_mesh()) m_have_mesh_channel=true;
                                        }
                                        else//continue with channel mesh calculations
                                        {
                                            if(m_have_mesh_channel)
                                            {
                                                if(!m_have_mesh_channel_normals) calc_mesh_channel_normals();
                                                else
                                                {
                                                    if(!m_have_real_mesh_channel_normals) calc_real_mesh_channel_normals();
                                                    else//done with channel mesh
                                                    {
                                                        if(calc_channel_mesh_vol()) //
                                                        {
                                                            HANDLE hConsole;
                                                            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                                                            SetConsoleTextAttribute(hConsole, 11); //set text color
                                                            cout<<"\nMesh Volume of Protein Channel is "<<m_channel_mesh_vol<<" A^3\n";
                                                            SetConsoleTextAttribute(hConsole, 15); //set text color
                                                        }
                                                        m_calc_done=true;//done with calculations
                                                        print_controls();
                                                        SetForegroundWindow(m_hWnd_gui); //make GUI window active
                                                    }
                                                }
                                            }
                                            else//without channel mesh
                                            {
                                                m_calc_done=true;//done with calculations
                                                print_controls();
                                                SetForegroundWindow(m_hWnd_gui); //make GUI window active
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(m_black_background)
    {
        glClearColor(0,0,0,1);
        GLfloat fogColor[4]= {0.01f, 0.01f, 0.01f, 0.5f};//Fog Color
        glFogfv(GL_FOG_COLOR, fogColor);//Set Fog Color
    }
    else
    {
        glClearColor(1,1,1,1);
        GLfloat fogColor[4]= {0.99f, 0.99f, 0.99f, 0.5f};//Fog Color
        glFogfv(GL_FOG_COLOR, fogColor);//Set Fog Color
    }

    glPushMatrix();
    //Drawing
    if(m_anaglyph) //for Red/Green 3D
    {
        //Stereo view
        stereo_cam cam(m_convergence,     // Convergence
                       m_eye_separation,       // Eye Separation
                       (float)m_window_length/(float)m_window_height,     // Aspect Ratio
                       45.0f,       // FOV along Y in degrees
                       0.10f,       // Near Clipping Distance
                       500.0f);   // Far Clipping Distance
        if(m_switch_eyes) cam.ApplyLeftFrustum();
        else cam.ApplyRightFrustum();

        set_camera_pos();

        glColorMask(true, false, false, false); //red layer
        draw_protein(); //draw protein
        glPopMatrix();glPushMatrix();
        glClear(GL_DEPTH_BUFFER_BIT) ;
        if(m_switch_eyes) cam.ApplyRightFrustum();
        else cam.ApplyLeftFrustum();

        set_camera_pos();

        //glTranslatef(0,0,1);
        glColorMask(false, true, true, false); //cyan layer
        draw_protein(); //Draw again
        glColorMask(true, true, true, true); //restore color
    }
    else if(m_interlaced) //for polarized 3D
    {
        //Fixing strips on stencil buffer
        GLint gliStencilBits;
        glGetIntegerv(GL_STENCIL_BITS,&gliStencilBits);
        int gliWindowWidth=m_window_length;
        int gliWindowHeight=m_window_height;
        //setting screen-corresponding geometry
        glViewport(0,0,gliWindowWidth,gliWindowHeight);
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity();
        gluOrtho2D(0.0,gliWindowWidth-1,0.0,gliWindowHeight-1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        //clearing and configuring stencil drawing
        glDrawBuffer(GL_BACK);
        glEnable(GL_STENCIL_TEST);
        glClearStencil(0);
        glClear(GL_STENCIL_BUFFER_BIT);
        glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE); // colorbuffer is copied to stencil
        glDisable(GL_DEPTH_TEST);
        glStencilFunc(GL_ALWAYS,1,1); // to avoid interaction with stencil content

        // drawing stencil pattern
        glColor4f(1,1,1,0);	// alfa is 0 not to interfere with alpha tests
        glLineWidth(1);
        glBegin(GL_LINES);
        for (int gliY=0; gliY<gliWindowHeight; gliY+=2)
        {
            glVertex2f(0,gliY);
            glVertex2f(gliWindowWidth,gliY);
        }
        glEnd();

        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP); // disabling changes in stencil buffer
        glStencilFunc(GL_NOTEQUAL,1,1); // draws if stencil is not 1
        glDrawBuffer(GL_BACK);
        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);

        //Stereo view
        stereo_cam cam(m_convergence,                                // Convergence
                       m_eye_separation,                             // Eye Separation
                       (float)m_window_length/(float)m_window_height,// Aspect Ratio
                       45.0f,                                        // FOV along Y in degrees
                       0.10f,                                        // Near Clipping Distance
                       500.0f);                                      // Far Clipping Distance
        //adjust position for one eye
        if(m_switch_eyes) cam.ApplyLeftFrustum();
        else cam.ApplyRightFrustum();

        set_camera_pos(); //place camera
        draw_protein(); //draw protein
        glPopMatrix(); //reset viewing matrix

        //Start with the other eye
        glPushMatrix();
        if(m_switch_eyes) cam.ApplyRightFrustum();
        else cam.ApplyLeftFrustum();

        set_camera_pos();
        glStencilFunc(GL_EQUAL,1,1);//only allow drawing on stencil buffer
        draw_protein(); //Draw again
        glDisable(GL_STENCIL_TEST);
    }
    else//not 3D
    {
        set_camera_pos();
        draw_protein();
    }
    glPopMatrix();
}


//Drawing

void protein::draw_protein(void)
{
    if(m_show_atoms)              draw_atoms();
    if(m_show_chains)             draw_alpha_carbons();
    if(m_show_channel_ruler)      draw_channel_edges();
    if(m_show_grid_protein)       draw_grid_all();
    if(m_show_grid_void)          draw_grid_void_volume();
    if(m_show_grid_void_exposure) draw_grid_void_color();
    if(m_show_other_atoms)        draw_other_atoms();
    //draw_grid_borders(); //debug
    if(m_show_channel_dist)       draw_convex_channel();
    if(m_have_mesh_atom&&m_show_mesh_atom) draw_atom_mesh();
    if(m_have_mesh_channel&&m_show_concave_channel) draw_concave_channel_mesh();
    if(m_show_labels)             draw_labels();
}

void protein::draw_atoms(void)
{
    //Icosahed shape
    float vertices[] = {
        0.276f,  0.851f,  0.447f, 0.894f,  0.000f,  0.447f, 0.000f,  0.000f,  1.000f,
        -0.724f,  0.526f,  0.447f, 0.276f,  0.851f,  0.447f, 0.000f,  0.000f,  1.000f,
        -0.724f, -0.526f,  0.447f, -0.724f,  0.526f,  0.447f, 0.000f,  0.000f,  1.000f,
        0.276f, -0.851f,  0.447f, -0.724f, -0.526f,  0.447f, 0.000f,  0.000f,  1.000f,
        0.894f,  0.000f,  0.447f, 0.276f, -0.851f,  0.447f, 0.000f,  0.000f,  1.000f,
        0.000f,  0.000f, -1.000f, 0.724f,  0.526f, -0.447f,  -0.276f,  0.851f, -0.447f,
        0.000f,  0.000f, -1.000f, -0.276f,  0.851f, -0.447f,  -0.894f,  0.000f, -0.447f,
        0.000f,  0.000f, -1.000f, -0.894f,  0.000f, -0.447f,  -0.276f, -0.851f, -0.447f,
        0.000f,  0.000f, -1.000f, -0.276f, -0.851f, -0.447f,  0.724f, -0.526f, -0.447f,
        0.000f,  0.000f, -1.000f, 0.724f, -0.526f, -0.447f, 0.724f,  0.526f, -0.447f,
        0.894f,  0.000f,  0.447f, 0.276f,  0.851f,  0.447f, 0.724f,  0.526f, -0.447f,
        0.276f,  0.851f,  0.447f, -0.724f,  0.526f,  0.447f, -0.276f,  0.851f, -0.447f,
        -0.724f,  0.526f,  0.447f, -0.724f, -0.526f,  0.447f, -0.894f,  0.000f, -0.447f,
        -0.724f, -0.526f,  0.447f, 0.276f, -0.851f,  0.447f, -0.276f, -0.851f, -0.447f,
        0.276f, -0.851f,  0.447f, 0.894f,  0.000f,  0.447f, 0.724f, -0.526f, -0.447f,
        0.276f,  0.851f,  0.447f,  -0.276f,  0.851f, -0.447f, 0.724f,  0.526f, -0.447f,
        -0.724f,  0.526f,  0.447f,  -0.894f,  0.000f, -0.447f, -0.276f,  0.851f, -0.447f,
        -0.724f, -0.526f,  0.447f,  -0.276f, -0.851f, -0.447f, -0.894f,  0.000f, -0.447f,
        0.276f, -0.851f,  0.447f, 0.724f, -0.526f, -0.447f, -0.276f, -0.851f, -0.447f,
        0.894f,  0.000f,  0.447f, 0.724f,  0.526f, -0.447f, 0.724f, -0.526f, -0.447f };
    float normals[] = {
        -0.275974,-0.850921,-0.446958,-0.894427,-0,-0.447214,-0,-0,-1,0.723761,-0.525826,
        -0.446852,-0.275974,-0.850921,-0.446958,-0,-0,-1,0.723761,0.525826,-0.446852,0.723761,
        -0.525826,-0.446852,-0,-0,-1,-0.275974,0.850921,-0.446958,0.723761,0.525826,-0.446852,
        -0,-0,-1,-0.894427,-0,-0.447214,-0.275974,0.850921,-0.446958,-0,-0,-1,-0,-0,1,-0.723761,
        -0.525826,0.446852,0.275974,-0.850921,0.446958,-0,-0,1,0.275974,-0.850921,0.446958,
        0.894427,-0,0.447214,-0,-0,1,0.894427,-0,0.447214,0.275974,0.850921,0.446958,-0,-0,
        1,0.275974,0.850921,0.446958,-0.723761,0.525826,0.446852,-0,-0,1,-0.723761,0.525826,
        0.446852,-0.723761,-0.525826,0.446852,-0.894427,-0,-0.447214,-0.275974,-0.850921,
        -0.446958,-0.723761,-0.525826,0.446852,-0.275974,-0.850921,-0.446958,0.723761,-0.525826,
        -0.446852,0.275974,-0.850921,0.446958,0.723761,-0.525826,-0.446852,0.723761,0.525826,
        -0.446852,0.894427,-0,0.447214,0.723761,0.525826,-0.446852,-0.275974,0.850921,-0.446958,
        0.275974,0.850921,0.446958,-0.275974,0.850921,-0.446958,-0.894427,-0,-0.447214,-0.723761,
        0.525826,0.446852,-0.275974,-0.850921,-0.446958,0.275974,-0.850921,0.446958,-0.723761,
        -0.525826,0.446852,0.723761,-0.525826,-0.446852,0.894427,-0,0.447214,0.275974,-0.850921,
        0.446958,0.723761,0.525826,-0.446852,0.275974,0.850921,0.446958,0.894427,-0,0.447214,
        -0.275974,0.850921,-0.446958,-0.723761,0.525826,0.446852,0.275974,0.850921,0.446958,
        -0.894427,-0,-0.447214,-0.723761,-0.525826,0.446852,-0.723761,0.525826,0.446852};

    //color
    GLfloat LightAmbientGray[]=		{ 0.3f, 0.3f, 0.4f, 0.9f };
    GLfloat LightDiffuseGray[]=		{ 0.5f, 0.5f, 0.5f, 0.9f };
    GLfloat LightSpecularGray[]=	{ 0.2f, 0.3f, 0.2f, 0.9f };

    GLfloat LightAmbientWhite[]=	{ 0.6f, 0.6f, 0.6f, 0.9f };
    GLfloat LightDiffuseWhite[]=	{ 0.8f, 0.8f, 0.8f, 0.9f };
    GLfloat LightSpecularWhite[]=	{ 0.9f, 0.9f, 0.9f, 0.9f };

    GLfloat LightAmbientRed[]=		{ 0.5f, 0.2f, 0.2f, 0.9f };
    GLfloat LightDiffuseRed[]=		{ 0.7f, 0.3f, 0.3f, 0.9f };
    GLfloat LightSpecularRed[]=	    { 0.3f, 0.2f, 0.2f, 0.9f };

    GLfloat LightAmbientBlue[]=		{ 0.3f, 0.3f, 0.5f, 0.9f };
    GLfloat LightDiffuseBlue[]=		{ 0.3f, 0.3f, 0.8f, 0.9f };
    GLfloat LightSpecularBlue[]=	{ 0.2f, 0.2f, 0.3f, 0.9f };

    GLfloat LightAmbientYellow[]=	{ 0.5f, 0.5f, 0.2f, 0.9f };
    GLfloat LightDiffuseYellow[]=	{ 0.5f, 0.5f, 0.1f, 0.9f };
    GLfloat LightSpecularYellow[]=	{ 0.3f, 0.3f, 0.2f, 0.9f };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);

    //drawing
    for (int i=0;i<(int)m_atoms.size();i++)
    {
        glPushMatrix();
        //set color
        switch(m_atoms[i].type)
        {
            case 'C':
            {
                glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbientGray);
                glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuseGray);
                glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecularGray);
            }break;
            case 'O':
            {
                glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbientRed);
                glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuseRed);
                glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecularRed);
            }break;
            case 'N':
            {
                glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbientBlue);
                glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuseBlue);
                glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecularBlue);
            }break;
            case 'S':
            {
                glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbientYellow);
                glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuseYellow);
                glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecularYellow);
            }break;
            case 'H':
            {
                glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbientWhite);
                glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuseWhite);
                glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecularWhite);
            }break;
        }

        glTranslatef(m_atoms[i].x,m_atoms[i].y,m_atoms[i].z);
        glScalef(m_atoms[i].radius,m_atoms[i].radius,m_atoms[i].radius);
        glDrawArrays(GL_TRIANGLES, 0, 60);
        glPopMatrix();
    }

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
}

void protein::draw_grid_all(void)
{
    if(m_pGrid==0) return;//have no grid yet

    glEnable(GL_DEPTH_TEST);
    glPointSize(3.0); //change size of a point
    glPushMatrix();
    glTranslatef(m_grid_pos.x,m_grid_pos.y,m_grid_pos.z);

    if(m_grid_updated) //draw grid if updated
    {
        glBegin(GL_POINTS);
        glColor4f(0.3,0.8,0.1,1);
        for(int x=0;x<m_grid_size[0];x++)
        {
            for(int y=0;y<m_grid_size[1];y++)
            {
                for(int z=0;z<m_grid_size[2];z++)
                {
                    if(m_pGrid[x][y][z])
                    {
                        glVertex3f((float)x*m_grid_gap_size,
                                   (float)y*m_grid_gap_size,
                                   (float)z*m_grid_gap_size);
                    }
                }
            }
        }
        glEnd();
    }
    else//draw progress in red
    {
        glBegin(GL_POINTS);
        glColor4f(1,0,0,1);
        for(int x=0;x<m_grid_size[0];x++)
        {
            for(int y=0;y<m_grid_size[1];y++)
            {
                for(int z=0;z<m_grid_z_progress;z++)
                {
                    if(m_pGrid[x][y][z])
                    {
                        glVertex3f((float)x*m_grid_gap_size,
                                   (float)y*m_grid_gap_size,
                                   (float)z*m_grid_gap_size);
                    }
                }
            }
        }
        glEnd();
    }
    glDisable(GL_DEPTH_TEST);
    glPopMatrix();
}

void protein::draw_grid_borders(void)
{
    glPushMatrix();
    glTranslatef(m_grid_pos.x,m_grid_pos.y,m_grid_pos.z);
    //axis
    glBegin(GL_LINES);
    glColor4f(1,0,0,1);
    glVertex3f(4.3,-0.8,-0.8); glVertex3f(-0.8,-0.8,-0.8);
    glColor4f(0,1,0,1);
    glVertex3f(-0.8,4.3,-0.8); glVertex3f(-0.8,-0.8,-0.8);
    glColor4f(0,0,1,1);
    glVertex3f(-0.8,-0.8,4.3); glVertex3f(-0.8,-0.8,-0.8);
    glEnd();
    glPointSize(1.0);
    //borders
    glBegin(GL_LINES);
    glColor4f(1,0,1,1);
    glVertex3f(-0.5,-0.5,-0.5); glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,-0.5);
    glVertex3f(-0.5,-0.5,-0.5); glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,-0.5);
    glVertex3f(-0.5,-0.5,-0.5); glVertex3f(-0.5,-0.5,m_grid_size[2]*m_grid_gap_size);

    glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,-0.5); glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,-0.5);
    glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,-0.5); glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,m_grid_size[2]*m_grid_gap_size);
    glVertex3f(-0.5,-0.5,m_grid_size[2]*m_grid_gap_size); glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);

    glVertex3f(-0.5,-0.5,m_grid_size[2]*m_grid_gap_size); glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,m_grid_size[2]*m_grid_gap_size);
    glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,-0.5); glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,-0.5);
    glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,-0.5); glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);

    glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,-0.5); glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);
    glVertex3f(-0.5,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size); glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);
    glVertex3f(m_grid_size[0]*m_grid_gap_size,-0.5,m_grid_size[2]*m_grid_gap_size); glVertex3f(m_grid_size[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);
    glEnd();
    glPopMatrix();
}

void protein::draw_grid_void_volume(void)
{
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    glTranslatef(m_grid_pos.x,m_grid_pos.y,m_grid_pos.z);

    if(m_void_vol_set) //draw grid if updated
    {
        //draw void blocks
        glColor3f(0.3,0.8,0.8);
        glPointSize(3);
        glBegin(GL_POINTS);
        for(int index=0;index<(int)m_void_blocks.size();index++)
        {
            glVertex3f((float)m_void_blocks[index].block_pos[0]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[1]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[2]*m_grid_gap_size);
        }
        glEnd();
    }
    else if(m_marking_mode)
    {
        //draw blocks
        //glDepthMask(GL_FALSE);
        for(int index=0;index<(int)m_void_blocks.size();index++)
        {

            if(m_void_blocks[index].block_pos[0]==m_marked_block[0]&&
               m_void_blocks[index].block_pos[1]==m_marked_block[1]&&
               m_void_blocks[index].block_pos[2]==m_marked_block[2])
            {
                glPointSize(10);
                glColor3f(1,0,0);
            }
            else if(m_void_blocks[index].block_pos[0]==m_marked_block[0]||
                    m_void_blocks[index].block_pos[1]==m_marked_block[1]||
                    m_void_blocks[index].block_pos[2]==m_marked_block[2])
            {
                glPointSize(6);
                glColor3f(0.5,0.5,1);
            }
            else
            {
                glPointSize(5);
                glColor3f(0,0,1);
            }
            glBegin(GL_POINTS);
            glVertex3f((float)m_void_blocks[index].block_pos[0]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[1]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[2]*m_grid_gap_size);
            glEnd();
        }

        //glDepthMask(GL_TRUE);
        //draw crosshair
        glLineWidth(2);
        glBegin(GL_LINES);
        //x
        glColor3f(0.2,0.0,0.0);
        glVertex3f(0,                                       (float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glColor3f(0.9,0.1,0.1);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glColor3f(1.0,0.7,0.7);
        glVertex3f(m_grid_size[0]*m_grid_gap_size,          (float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        //y
        glColor3f(0.0,0.2,0.0);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,0                                       ,(float)m_marked_block[2]*m_grid_gap_size);
        glColor3f(0.1,0.9,0.1);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glColor3f(0.7,1.0,0.7);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,m_grid_size[1]*m_grid_gap_size,          (float)m_marked_block[2]*m_grid_gap_size);
        //z
        glColor3f(0.1,0.1,0.3);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,0                                       );
        glColor3f(0.2,0.2,0.9);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,(float)m_marked_block[2]*m_grid_gap_size);
        glColor3f(0.7,0.7,1.0);
        glVertex3f((float)m_marked_block[0]*m_grid_gap_size,(float)m_marked_block[1]*m_grid_gap_size,m_grid_size[2]*m_grid_gap_size);
        glEnd();

        if(m_channel_is_leaking) //for stopping leakage
        {
            //draw marked blocks
            glPointSize(12);
            glColor3f(0,1,0);
            glBegin(GL_POINTS);
            for(int index=0;index<(int)m_void_blocks_marked.size();index++)
            {

                glVertex3f((float)m_void_blocks_marked[index].block_pos[0]*m_grid_gap_size,
                           (float)m_void_blocks_marked[index].block_pos[1]*m_grid_gap_size,
                           (float)m_void_blocks_marked[index].block_pos[2]*m_grid_gap_size);
            }
            glEnd();
        }
        else //for marking channel
        {
            glPointSize(12);
            if(m_channel_start_set)
            {
                //draw start block
                glColor3f(0.9,0.9,0.3);
                glBegin(GL_POINTS);
                glVertex3f(m_channel_start_pos.x*m_grid_gap_size,
                           m_channel_start_pos.y*m_grid_gap_size,
                           m_channel_start_pos.z*m_grid_gap_size);
                glEnd();
            }
            if(m_channel_end_set)
            {
                //draw end block
                glColor3f(0.9,0.9,0.1);
                glBegin(GL_POINTS);
                glVertex3f(m_channel_end_pos.x*m_grid_gap_size,
                           m_channel_end_pos.y*m_grid_gap_size,
                           m_channel_end_pos.z*m_grid_gap_size);
                glEnd();
            }
            if(m_channel_center_set)
            {
                //draw center block
                glColor3f(0.9,0.7,0.1);
                glBegin(GL_POINTS);
                glVertex3f(m_channel_center_pos.x*m_grid_gap_size,
                           m_channel_center_pos.y*m_grid_gap_size,
                           m_channel_center_pos.z*m_grid_gap_size);
                glEnd();
            }
            if(m_channel_start_set&&m_channel_end_set)
            {
                //draw channel axis
                glLineWidth(5);
                glColor3f(0,1,1);
                glBegin(GL_LINES);
                glVertex3f(m_channel_start_pos.x*m_grid_gap_size,
                           m_channel_start_pos.y*m_grid_gap_size,
                           m_channel_start_pos.z*m_grid_gap_size);
                glVertex3f(m_channel_end_pos.x*m_grid_gap_size,
                           m_channel_end_pos.y*m_grid_gap_size,
                           m_channel_end_pos.z*m_grid_gap_size);
                glEnd();
            }
        }
    }
    else//draw progress
    {
        //draw void blocks
        glPointSize(3);
        glBegin(GL_POINTS);
        glColor4f(1,0,1,1);
        for(int index=0;index<(int)m_void_blocks.size();index++) //draw stored blocks
        {
            glVertex3f((float)m_void_blocks[index].block_pos[0]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[1]*m_grid_gap_size,
                       (float)m_void_blocks[index].block_pos[2]*m_grid_gap_size);
        }
        glColor4f(1,1,1,1);
        for(int index=0;index<(int)m_void_blocks_new.size();index++) //and the new ones
        {
            glVertex3f((float)m_void_blocks_new[index].block_pos[0]*m_grid_gap_size,
                       (float)m_void_blocks_new[index].block_pos[1]*m_grid_gap_size,
                       (float)m_void_blocks_new[index].block_pos[2]*m_grid_gap_size);
        }
        glEnd();
    }
    //glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glPopMatrix();
}

void protein::draw_grid_void_color(void)
{
    if(m_pGrid_color&&m_have_grid_void_color)//cant draw if grid have no grid
    {
        glPointSize(5);
        glPushMatrix();
        glTranslatef(m_grid_pos.x,m_grid_pos.y,m_grid_pos.z);
        glEnable(GL_DEPTH_TEST);
        glBegin(GL_POINTS);
        for(int x=0;x<m_grid_size[0];x++)
        for(int y=0;y<m_grid_size[1];y++)
        for(int z=0;z<m_grid_size[2];z++)
        {
            if(m_show_grid_void_exposure_level<m_pGrid_color[x][y][z]) continue; //show specified exposure
            switch(m_pGrid_color[x][y][z])//set color depending on exposure
            {
                case 0: glColor3f(1.0,0.0,0.0);break;
                case 1: glColor3f(0.9,0.4,0.4);break;
                case 2: glColor3f(0.7,0.2,0.2);break;
                case 3: glColor3f(0.5,0.1,0.1);break;
                case 4: glColor3f(0.3,0.2,0.1);break;
                case 5: glColor3f(0.1,0.3,0.0);break;
                case 6: glColor3f(0.0,0.5,0.0);break;
                case 7: continue;//skip if part of protein
                default:glColor3f(1.0,1.0,1.0);break;
            }
            glVertex3f((float)x*m_grid_gap_size,
                       (float)y*m_grid_gap_size,
                       (float)z*m_grid_gap_size);
        }
        glEnd();
        glDisable(GL_DEPTH_TEST);
        glPopMatrix();
    }
}

void protein::draw_alpha_carbons(void)
{
    glLineWidth(3);
    glEnable(GL_DEPTH_TEST);
    glPushMatrix();
    for(int chain=0;chain<(int)m_chains.size();chain++)
    {
        switch(chain)//set color
        {
            case  0: glColor3f(0.9,0.3,0.1);break;
            case  1: glColor3f(0.3,0.9,0.1);break;
            case  2: glColor3f(0.3,0.1,0.9);break;
            case  3: glColor3f(0.9,0.3,0.9);break;
            case  4: glColor3f(0.3,0.9,0.3);break;
            case  5: glColor3f(0.0,0.7,0.7);break;
            case  6: glColor3f(0.4,0.4,0.6);break;
            case  7: glColor3f(0.8,0.8,0.0);break;
            case  8: glColor3f(0.8,0.4,0.4);break;
            case  9: glColor3f(0.0,0.9,0.6);break;
            default: glColor3f(0.9,0.9,0.9);break;
        }
        glBegin(GL_LINE_STRIP);
        for(int i=0;i<(int)m_chains[chain].size();i++)
        {
            glVertex3f(m_chains[chain][i].x,m_chains[chain][i].y,m_chains[chain][i].z);
        }
        glEnd();
    }
    glPopMatrix();
    glDisable(GL_DEPTH_TEST);
    glLineWidth(1);
}

void protein::draw_channel_edges(void)
{
    if(m_channel_edges_up.empty()) return;
    int roll=m_axis_roll;
    while(roll>=360) roll-=360;
    while(roll<0) roll+=360;
    bool col=false;
    if(roll>=180)
    {
        roll-=180;
    }
    glEnable(GL_DEPTH_TEST);
    if(m_channel_edge_set) //calculations done
    {
        glColor3f(0,1,0);
        glBegin(GL_LINES);
        for(int point=0;point<(int)m_channel_edges_up[roll].size();point++)
        {
            glVertex3f(m_channel_edges_up[int(roll)][point].x,
                       m_channel_edges_up[int(roll)][point].y,
                       m_channel_edges_up[int(roll)][point].z);
            glVertex3f(m_channel_edges_up[int(roll+180)][point].x,
                       m_channel_edges_up[int(roll+180)][point].y,
                       m_channel_edges_up[int(roll+180)][point].z);
        }
        glEnd();
    }
    else if(m_axis_roll_progress>0) //draw roll that just been calculated
    {
        //draw both up and down
        roll=m_axis_roll_progress-1;
        glColor3f(1,0.5,0);
        glBegin(GL_LINES);
        for(int point=0;point<(int)m_channel_edges_up[roll].size();point++)
        {
            glVertex3f(m_channel_edges_up[int(roll)][point].x,
                       m_channel_edges_up[int(roll)][point].y,
                       m_channel_edges_up[int(roll)][point].z);
            glVertex3f(m_channel_edges_down[int(roll)][point].x,
                       m_channel_edges_down[int(roll)][point].y,
                       m_channel_edges_down[int(roll)][point].z);
        }
        glEnd();
    }
    glDisable(GL_DEPTH_TEST);
}

void protein::draw_other_atoms(void)
{
    glEnable(GL_DEPTH_TEST);
    glPointSize(3.0);
    glPushMatrix();
    glBegin(GL_POINTS);
    glColor3f(1,1,1);
    for(int atom=0;atom<(int)m_atoms_other.size();atom++)
    {
        glVertex3f(m_atoms_other[atom].x,m_atoms_other[atom].y,m_atoms_other[atom].z);
    }
    glEnd();
    glPopMatrix();
    glPointSize(1.0);
    glDisable(GL_DEPTH_TEST);
}

void protein::draw_convex_channel(void)
{
    if(m_channel_edges_up.empty()) return;//no drawing if no data
    glEnable(GL_DEPTH_TEST);
    glColor3f(0.7,0.8,0.2);//loading color
    //draw line along the roll axis
    int number_of_segments=(int)m_channel_edges_up[0].size();//ok if all rolls have same number of segments
    for(int segment=0;segment<number_of_segments;segment++)
    {
        glBegin(GL_LINE_STRIP);
        for(int roll=0;roll<(int)m_channel_edges_up.size();roll++)
        {
            if(m_have_convex_channel_lengths)
            {
                glColor3f(m_channel_edges_color[roll][segment].x,
                          m_channel_edges_color[roll][segment].y,
                          m_channel_edges_color[roll][segment].z);
            }
            glVertex3f(m_channel_edges_up[roll][segment].x,
                       m_channel_edges_up[roll][segment].y,
                       m_channel_edges_up[roll][segment].z);
        }
        glEnd();
    }
    //draw line along the segment axis
    for(int roll=0;roll<(int)m_channel_edges_up.size();roll++)
    {
        glBegin(GL_LINE_STRIP);
        for(int segment=0;segment<(int)m_channel_edges_up[roll].size();segment++)
        {
            if(m_have_convex_channel_lengths)
            {
                glColor3f(m_channel_edges_color[roll][segment].x,
                          m_channel_edges_color[roll][segment].y,
                          m_channel_edges_color[roll][segment].z);
            }
            glVertex3f(m_channel_edges_up[roll][segment].x,
                       m_channel_edges_up[roll][segment].y,
                       m_channel_edges_up[roll][segment].z);
        }
        glEnd();
    }
    glDisable(GL_DEPTH_TEST);
}

void protein::draw_atom_mesh(void)
{
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    if(m_mesh_transparency) glEnable(GL_BLEND);
    //set color
    GLfloat LightAmbient[]=		{ 0.3f, 0.3f, 0.4f, 0.4f };
    GLfloat LightDiffuse[]=		{ 0.7f, 0.7f, 0.8f, 0.4f };
    GLfloat LightSpecular[]=	{ 0.2f, 0.3f, 0.2f, 0.3f };
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbient);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuse);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecular);
    glEnable(GL_DEPTH_TEST);
    if(m_have_real_mesh_normals)//draw with real normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use average normal
        for(int i=0;i<(int)m_mesh_atom_triangles.size();i++)
        {
            glNormal3f(m_mesh_atom_triangles[i].first_avg_normal.x,
                       m_mesh_atom_triangles[i].first_avg_normal.y,
                       m_mesh_atom_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].z,1);
            glNormal3f(m_mesh_atom_triangles[i].second_avg_normal.x,
                       m_mesh_atom_triangles[i].second_avg_normal.y,
                       m_mesh_atom_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].z,1);
            glNormal3f(m_mesh_atom_triangles[i].third_avg_normal.x,
                       m_mesh_atom_triangles[i].third_avg_normal.y,
                       m_mesh_atom_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else if(m_have_mesh_normals)//draw with common normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_vertex_progress;i++)
        {
            glNormal3f(m_mesh_atom_triangles[i].first_avg_normal.x,
                       m_mesh_atom_triangles[i].first_avg_normal.y,
                       m_mesh_atom_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].z,1);
            glNormal3f(m_mesh_atom_triangles[i].second_avg_normal.x,
                       m_mesh_atom_triangles[i].second_avg_normal.y,
                       m_mesh_atom_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].z,1);
            glNormal3f(m_mesh_atom_triangles[i].third_avg_normal.x,
                       m_mesh_atom_triangles[i].third_avg_normal.y,
                       m_mesh_atom_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].z,1);
        }
        for(int i=m_mesh_normal_vertex_progress;i<m_mesh_normal_progress;i++) //draw rest with common normals
        {
            glNormal3f(m_mesh_atom_triangles[i].first_normal.x,
                       m_mesh_atom_triangles[i].first_normal.y,
                       m_mesh_atom_triangles[i].first_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].z,1);
            glNormal3f(m_mesh_atom_triangles[i].second_normal.x,
                       m_mesh_atom_triangles[i].second_normal.y,
                       m_mesh_atom_triangles[i].second_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].z,1);
            glNormal3f(m_mesh_atom_triangles[i].third_normal.x,
                       m_mesh_atom_triangles[i].third_normal.y,
                       m_mesh_atom_triangles[i].third_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else
    {
        //draw triangles with normals calculated
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_progress;i++)
        {
            glNormal3f(m_mesh_atom_triangles[i].first_normal.x,
                       m_mesh_atom_triangles[i].first_normal.y,
                       m_mesh_atom_triangles[i].first_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].first].z,1);
            glNormal3f(m_mesh_atom_triangles[i].second_normal.x,
                       m_mesh_atom_triangles[i].second_normal.y,
                       m_mesh_atom_triangles[i].second_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].second].z,1);
            glNormal3f(m_mesh_atom_triangles[i].third_normal.x,
                       m_mesh_atom_triangles[i].third_normal.y,
                       m_mesh_atom_triangles[i].third_normal.z);
            glVertex4f(m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].x,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].y,
                       m_mesh_atom_vertices[m_mesh_atom_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    //glDisable(GL_POLYGON_SMOOTH);
}

void protein::draw_concave_channel_mesh(void)
{
    //glDepthMask(GL_FALSE);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glCullFace(GL_BACK);
    if(m_mesh_transparency)
    {
        glEnable(GL_BLEND);
        glEnable(GL_CULL_FACE);
        glDisable(GL_FOG);
        //set color
        GLfloat LightAmbient[]=		{ 0.3f, 0.0f, 0.0f, 0.5f };
        GLfloat LightDiffuse[]=		{ 0.9f, 0.0f, 0.0f, 0.5f };
        GLfloat LightSpecular[]=	{ 0.9f, 0.5f, 0.5f, 0.5f };
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbient);
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuse);
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecular);
    }
    else
    {
        if(m_channel_mesh_see_through)
        {
            glEnable(GL_CULL_FACE);
            glCullFace(GL_FRONT);
        }
        glEnable(GL_DEPTH_TEST);
        //set color
        GLfloat LightAmbient[]=		{ 0.5f, 0.1f, 0.1f, 1.0f };
        GLfloat LightDiffuse[]=		{ 0.9f, 0.3f, 0.3f, 1.0f };
        GLfloat LightSpecular[]=	{ 0.9f, 0.3f, 0.3f, 1.0f };
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,LightAmbient);
        glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,LightDiffuse);
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,LightSpecular);
    }
    //glEnable(GL_DEPTH_TEST);

    if(m_have_real_mesh_channel_normals)//draw with real normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use average normal
        for(int i=0;i<(int)m_mesh_channel_triangles.size();i++)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_avg_normal.x,
                       m_mesh_channel_triangles[i].first_avg_normal.y,
                       m_mesh_channel_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_avg_normal.x,
                       m_mesh_channel_triangles[i].second_avg_normal.y,
                       m_mesh_channel_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_avg_normal.x,
                       m_mesh_channel_triangles[i].third_avg_normal.y,
                       m_mesh_channel_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else if(m_have_mesh_channel_normals)//draw with common normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_vertex_progress;i++)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_avg_normal.x,
                       m_mesh_channel_triangles[i].first_avg_normal.y,
                       m_mesh_channel_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_avg_normal.x,
                       m_mesh_channel_triangles[i].second_avg_normal.y,
                       m_mesh_channel_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_avg_normal.x,
                       m_mesh_channel_triangles[i].third_avg_normal.y,
                       m_mesh_channel_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        for(int i=m_mesh_normal_vertex_progress;i<m_mesh_normal_progress;i++) //draw rest with common normals
        {
            glNormal3f(m_mesh_channel_triangles[i].first_normal.x,
                       m_mesh_channel_triangles[i].first_normal.y,
                       m_mesh_channel_triangles[i].first_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_normal.x,
                       m_mesh_channel_triangles[i].second_normal.y,
                       m_mesh_channel_triangles[i].second_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_normal.x,
                       m_mesh_channel_triangles[i].third_normal.y,
                       m_mesh_channel_triangles[i].third_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else
    {
        //draw triangles with normals calculated
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_progress;i++)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_normal.x,
                       m_mesh_channel_triangles[i].first_normal.y,
                       m_mesh_channel_triangles[i].first_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_normal.x,
                       m_mesh_channel_triangles[i].second_normal.y,
                       m_mesh_channel_triangles[i].second_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_normal.x,
                       m_mesh_channel_triangles[i].third_normal.y,
                       m_mesh_channel_triangles[i].third_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    if(!m_mesh_transparency) return;
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //BACK TO FRONT

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    if(m_mesh_transparency) glEnable(GL_BLEND);
    //glEnable(GL_DEPTH_TEST);

    if(m_have_real_mesh_channel_normals)//draw with real normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use average normal
        for(int i=(int)m_mesh_channel_triangles.size()-1;i>=0;i--)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_avg_normal.x,
                       m_mesh_channel_triangles[i].first_avg_normal.y,
                       m_mesh_channel_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_avg_normal.x,
                       m_mesh_channel_triangles[i].second_avg_normal.y,
                       m_mesh_channel_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_avg_normal.x,
                       m_mesh_channel_triangles[i].third_avg_normal.y,
                       m_mesh_channel_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else if(m_have_mesh_channel_normals)//draw with common normals
    {
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_vertex_progress;i++)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_avg_normal.x,
                       m_mesh_channel_triangles[i].first_avg_normal.y,
                       m_mesh_channel_triangles[i].first_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_avg_normal.x,
                       m_mesh_channel_triangles[i].second_avg_normal.y,
                       m_mesh_channel_triangles[i].second_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_avg_normal.x,
                       m_mesh_channel_triangles[i].third_avg_normal.y,
                       m_mesh_channel_triangles[i].third_avg_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        for(int i=m_mesh_normal_vertex_progress;i<m_mesh_normal_progress;i++) //draw rest with common normals
        {
            glNormal3f(m_mesh_channel_triangles[i].first_normal.x,
                       m_mesh_channel_triangles[i].first_normal.y,
                       m_mesh_channel_triangles[i].first_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_normal.x,
                       m_mesh_channel_triangles[i].second_normal.y,
                       m_mesh_channel_triangles[i].second_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_normal.x,
                       m_mesh_channel_triangles[i].third_normal.y,
                       m_mesh_channel_triangles[i].third_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    else
    {
        //draw triangles with normals calculated
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        //use common normals
        for(int i=0;i<m_mesh_normal_progress;i++)
        {
            glNormal3f(m_mesh_channel_triangles[i].first_normal.x,
                       m_mesh_channel_triangles[i].first_normal.y,
                       m_mesh_channel_triangles[i].first_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].first].z,1);
            glNormal3f(m_mesh_channel_triangles[i].second_normal.x,
                       m_mesh_channel_triangles[i].second_normal.y,
                       m_mesh_channel_triangles[i].second_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].second].z,1);
            glNormal3f(m_mesh_channel_triangles[i].third_normal.x,
                       m_mesh_channel_triangles[i].third_normal.y,
                       m_mesh_channel_triangles[i].third_normal.z);
            glVertex4f(m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].x,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].y,
                       m_mesh_channel_vertices[m_mesh_channel_triangles[i].third].z,1);
        }
        glEnd();
        glDisable(GL_LIGHTING);
    }
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glEnable(GL_FOG);

    //glDepthMask(GL_TRUE);

}

void protein::draw_labels(void)
{
    glColor3f(1,1,1);
    glDisable(GL_FOG);
    string label;
    for(int chain=0;chain<(int)m_chains.size();chain++)
    {
        for(int aatom=0;aatom<(int)m_chains[chain].size();aatom++)
        {
            float pos[]={m_chains[chain][aatom].x,m_chains[chain][aatom].y,m_chains[chain][aatom].z};
            //switch(m_chains[chain][aatom].type)
            label=m_chains[chain][aatom].residue;
            char buffer[10];
            itoa(aatom+1,buffer,10);
            label.append(buffer);
            switch(m_view_mode)
            {
                case view_arcball:
                {
                    float matrix[16];
                    m_arcball.get_inverse_matrix(matrix);
                    Text_box.draw_text(0.5,pos,label,matrix);
                }break;
                case view_locked_roll:
                {
                    vec channel_axis,rot_axis;
                    float rotA;
                    if(m_channel_defined)
                    {
                        //Normalize channel vector
                        channel_axis=m_channel_end_pos-m_channel_start_pos;
                        channel_axis=channel_axis.unit();
                        //get angle between vectors (dot product between channel and z-axis)
                        rotA=acos(channel_axis*vec(0,0,1))/_piover180;
                        //get perpendicular axis
                        rot_axis=channel_axis^vec(0,0,1);
                        rot_axis=rot_axis.unit();
                    }
                    else//Z-axis in channel vector
                    {
                        channel_axis=vec(0,0,1);
                        rotA=180;
                        rot_axis=vec(0,0,1);
                    }

                    Text_box.draw_text(0.5,pos,label,rot_axis,channel_axis,rotA,m_axis_roll);
                }break;
                case view_FPS:
                {
                    Text_box.draw_text(0.5,pos,label,m_rot_hori,m_rot_vert);
                }break;
            }
        }
    }
    glEnable(GL_FOG);
}

//Output

bool protein::output_diameters_along_axis(void)
{
    if(!m_want_dist_channel)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set text color
        cout<<"\nChannel have not been measured.\nReload protein if these values should be calculated.\n";
        SetConsoleTextAttribute(hConsole, 14); //set text color
        return false;
    }
    //create filename
    int roll=m_axis_roll;
    while(roll>360) roll-=360;
    while(roll<0) roll+=360;
    stringstream out;
    out << roll; //int to string conversion
    string sVal=out.str();
    string filename="OUTPUT\\Diameters_";
    filename.append(sVal);
    filename.append("_deg.txt");
    //create file
    ofstream file(filename.c_str());
    if(file==0) //error
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR writing file\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //write data
    file<<"Diameters of Protein Channel for Angle "<<roll;//Header
    file<<"\nIndex Diameter (Å)\n";
    if(roll>=180) roll-=180;
    for(int index=0;index<(int)m_channel_edges_length[roll].size();index++)
    {
        float diameter=m_channel_edges_length[roll][index]+m_channel_edges_length[roll+180][index];
        file<<index+1<<' ';
        file<<diameter<<endl;
    }
    file.close();
    string ofilename="Diameters_";
    ofilename.append(sVal);
    ofilename.append("_deg.txt");
    cout<<"File output: "<<ofilename<<endl;
    return true;
}

bool protein::output_diameters_average(void)
{
    if(!m_want_dist_channel)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set text color
        cout<<"\nChannel have not been measured.\nReload protein if these values should be calculated.\n";
        SetConsoleTextAttribute(hConsole, 14); //set text color
        return false;
    }
    //create file
    ofstream file("OUTPUT\\Diameters_average.txt");
    if(file==0) //error
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR writing file\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //write data
    file<<"Average Diameter of the Protein Channel Along the Channel";//Header
    file<<"\nIndex Diameter (Å)\n";
    for(int index=0;index<(int)m_channel_edges_length[0].size();index++)//assume that all angles have same number of measurments
    {
        //calc average radius
        float diameter_sum=0;
        for(int angle=0;angle<180;angle++)
        {
            float diameter=m_channel_edges_length[angle][index]+m_channel_edges_length[angle+180][index];
            diameter_sum+=diameter;
        }
        float diameter_average=diameter_sum/180;
        //print
        file<<index+1<<' ';
        file<<diameter_average<<endl;
    }
    file.close();
    cout<<"File output: Diameters_average.txt"<<endl;
    return true;
}

bool protein::output_diameters_min(void)
{
    if(!m_want_dist_channel)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set text color
        cout<<"\nChannel have not been measured.\nReload protein if these values should be calculated.\n";
        SetConsoleTextAttribute(hConsole, 14); //set text color
        return false;
    }
    //create file
    ofstream file("OUTPUT\\Diameters_minimum.txt");
    if(file==0) //error
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR writing file\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //write data
    file<<"Smallest Diameters of the Protein Channel Along the Channel";//Header
    file<<"\nIndex Diameter (Å)\n";
    for(int index=0;index<(int)m_channel_edges_length[0].size();index++)//asume all angles have same number of measurments
    {
        //find shortest distance
        float diameter_min=m_channel_edges_length[0][index]+m_channel_edges_length[180][index];
        for(int angle=0;angle<180;angle++)
        {
            float diameter=m_channel_edges_length[angle][index]+m_channel_edges_length[angle+180][index];
            if(diameter<diameter_min)
            {
                diameter_min=diameter;
            }
        }
        //print
        file<<index+1<<' ';
        file<<diameter_min<<endl;
    }
    file.close();
    cout<<"File output: Diameters_minimum.txt"<<endl;
    return true;
}

bool protein::output_radiuses_min(void)
{
    if(!m_want_dist_channel)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set text color
        cout<<"\nChannel have not been measured.\nReload protein if these values should be calculated.\n";
        SetConsoleTextAttribute(hConsole, 14); //set text color
        return false;
    }
    //create file
    ofstream file("OUTPUT\\Radiuses_minimum.txt");
    if(file==0) //error
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR writing file\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //write data
    file<<"Smallest Radiuses of the Protein Channel Along the Channel";//Header
    file<<"\nIndex Radius (Å)\n";
    for(int index=0;index<(int)m_channel_edges_length[0].size();index++)//asume all angles have same number of measurments
    {
        //find shortest distance
        float radius_min=m_channel_edges_length[0][index];
        for(int angle=0;angle<360;angle++)
        {
            if(m_channel_edges_length[angle][index]<radius_min)
            {
                radius_min=m_channel_edges_length[angle][index];
            }
        }
        //print
        file<<index+1<<' ';
        file<<radius_min<<endl;
    }
    file.close();
    cout<<"File output: Radiuses_minimum.txt"<<endl;
    return true;
}


//Misc

bool protein::move(bool keys[],bool mouse_buttons[],int mouse_pos[],float cycletime)
{
    if(keys[82]) //r
    {
        m_key_delay=200;
        //load a new protein
        cout<<string(20, '\n'); //clear screen
        SetForegroundWindow(m_hWnd_console); //make console window active
        return false;
    }

    m_key_delay-=cycletime; //update key delay timer

    //Output
    if (keys[80] && m_key_delay<=0 && m_calc_done) //p
    {
        output_diameters_along_axis();
        m_key_delay=1000;
    }
    if (keys[79] && m_key_delay<=0 && m_calc_done) //o
    {
        output_diameters_average();
        m_key_delay=1000;
    }
    if (keys[73] && m_key_delay<=0 && m_calc_done) //i
    {
        output_diameters_min();
        m_key_delay=1000;
    }
    if (keys[85] && m_key_delay<=0 && m_calc_done) //u
    {
        output_radiuses_min();
        m_key_delay=1000;
    }
    //view mode
    if (keys[9] && m_key_delay<=0) //Tab change viewmode
    {
        //cycle views
        if(m_view_mode==view_FPS) m_view_mode=view_locked_roll;
        else if(m_view_mode==view_locked_roll)
        {
            m_view_mode=view_arcball;
            init_arcball();
        }
        else
        {
            HANDLE hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
            SetConsoleTextAttribute(hConsole, 14); //set text color
            cout<<"\nFirst Person Camera Movment\nControls:\n"
                  "[W/S]   - Move Forward/Backward\n"
                  "[A/D]   - Move Left/Right\n"
                  "[Q/E]   - Move Up/Down\n"
                  "[Mouse] - Left-Click on the screen and drag to change viewing direction\n\n";
            SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
            //get position from arcball rotation
            m_eye.x=m_cam_pos.x;
            m_eye.y=m_cam_pos.y;
            m_eye.z=m_cam_pos.z;

            m_view_mode=view_FPS;
            m_prev_mouse_pos[0]=mouse_pos[0];
            m_prev_mouse_pos[1]=mouse_pos[1];
        }
        m_key_delay=200;
    }
    //Anaglyph 3D
    if(keys[113]&&m_key_delay<=0) //F2
    {
        m_anaglyph=!m_anaglyph;
        m_interlaced=false;
        m_key_delay=200;
    }
    //Interlaced 3D
    if(keys[114]&&m_key_delay<=0) //F3
    {
        m_interlaced=!m_interlaced;
        m_anaglyph=false;
        m_key_delay=200;
    }
    //Eye switch
    if(keys[115]&&m_key_delay<=0) //F4
    {
        m_switch_eyes=!m_switch_eyes;
        m_key_delay=200;
    }
    //3D controls
    if(m_anaglyph||m_interlaced)
    {
        if(keys[39]&&m_key_delay<=0) m_eye_separation+=0.05; //right
        if(keys[37]&&m_key_delay<=0) m_eye_separation-=0.05; //left
        if(keys[38]&&m_key_delay<=0) m_convergence+=0.1; //up
        if(keys[40]&&m_key_delay<=0) m_convergence-=0.1; //down
    }
    //Background color
    if(keys[66]&&m_key_delay<=0)//b
    {
        m_black_background=!m_black_background;
        Text_box.swap_background_color();
        m_key_delay=200;
    }
    //Toogle Transparency
    if(keys[84]&&m_key_delay<=0) //T
    {
        m_mesh_transparency=!m_mesh_transparency;
        m_key_delay=200;
    }
    //Show labels
    if(keys[76]&&m_key_delay<=0) //L
    {
        m_show_labels=!m_show_labels;
        m_key_delay=200;
    }
    //Define Channel Position
    if(keys[70]&&m_key_delay<=0) //F
    {
        m_axis_roll=0; //resets view
        m_axis_roll_progress=0;
        m_marked_block[0]=m_grid_size[0]/2;
        m_marked_block[1]=m_grid_size[1]/2;
        m_marked_block[2]=m_grid_size[2]/2;
        m_show_grid_void=true;
        m_channel_start_set=false;
        m_channel_end_set=false;
        m_channel_center_set=false;
        m_calc_done=false;
        m_want_defined_channel=true;
        m_channel_defined=false;
        m_void_vol_set=false;
        m_channel_edge_set=false;
        m_have_convex_channel_lengths=false;
        m_tried_to_get_channel_mesh=false;
        m_have_mesh_channel=false;
        m_have_mesh_channel_normals=false;
        m_have_real_mesh_channel_normals=false;
        m_axis_roll_progress=0;
        m_mesh_normal_progress=0;
        //clear data
        m_void_blocks.clear();
        m_void_blocks_exp.clear();
        m_void_blocks_new.clear();
        m_void_blocks_marked.clear();
        m_channel_edges_up.clear();
        m_channel_edges_normals.clear();
        m_channel_edges_color.clear();
        m_channel_edges_length.clear();
        m_mesh_channel_vertices.clear();
        m_mesh_channel_triangles.clear();
        m_vertex_triangle_map_channel.clear();
        //output
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 14); //set text color
        cout<<"You have selected to define the position and axis of the channel.\n"
              "Controls: X-Axis (RED), Y-Axis (GREEN), Z-Axis (BLUE)\n"
              "[Numpad 4] - Move marker in positive X direction\n"
              "[Numpad 1] - Move marker in negative X direction\n"
              "[Numpad 5] - Move marker in positive Y direction\n"
              "[Numpad 2] - Move marker in negative Y direction\n"
              "[Numpad 6] - Move marker in positive Z direction\n"
              "[Numpad 3] - Move marker in negative Z direction\n"
              "[Space]    - Mark selected grid block\n"
              "[Del]      - Remove selection\n"
              "[Enter]    - When marking is complete\n"
              "\nNow select the starting point of the channel\n";
        SetConsoleTextAttribute(hConsole, 15); //set text color
        m_key_delay=200;
    }
    //See-through channel mesh
    if(keys[71]&&m_key_delay<=0)//G
    {
        m_channel_mesh_see_through=!m_channel_mesh_see_through;
        m_key_delay=200;
    }
    //protein view
    if(!m_leakage_waiting_for_response)
    {
        if(keys[49] && m_key_delay<=0) {m_key_delay=200; m_show_atoms=!m_show_atoms;} //1
        if(keys[50] && m_key_delay<=0) {m_key_delay=200; m_show_chains=!m_show_chains;} //2
        if(keys[51] && m_key_delay<=0) {m_key_delay=200; m_show_grid_protein=!m_show_grid_protein;} //3
        if(keys[52] && m_key_delay<=0) {m_key_delay=200; m_show_grid_void=!m_show_grid_void;} //4
        if(keys[53] && m_key_delay<=0) {m_key_delay=200; m_show_channel_dist=!m_show_channel_dist;} //5
        if(keys[54] && m_key_delay<=0) {m_key_delay=200; m_show_channel_ruler=!m_show_channel_ruler;} //6
        if(keys[55] && m_key_delay<=0) {m_key_delay=200; m_show_mesh_atom=!m_show_mesh_atom;} //7
        if(keys[56] && m_key_delay<=0) {m_key_delay=200; m_show_concave_channel=!m_show_concave_channel;} //8
        if(keys[57] && m_key_delay<=0) {m_key_delay=200; m_show_grid_void_exposure=!m_show_grid_void_exposure;} //9
        if(keys[48] && m_key_delay<=0) {m_key_delay=200; m_show_other_atoms=!m_show_other_atoms;} //0
    }
    //Grid void exposure level
    if(keys[33]&&m_key_delay<=0&&m_show_grid_void_exposure_level<7) {m_show_grid_void_exposure_level++;m_key_delay=200;}
    if(keys[34]&&m_key_delay<=0&&m_show_grid_void_exposure_level>0) {m_show_grid_void_exposure_level--;m_key_delay=200;}

    //channel axis roll
    if(keys[107]) {m_axis_roll+=cycletime*0.1;} //+
    if(keys[109]) {m_axis_roll-=cycletime*0.1;} //-
    //reset viewing pos
    if(keys[8]) //backspace
    {
        m_eye.x=0; m_eye.y=0; m_eye.z=50;
        m_rot_hori=270; m_rot_vert=90;
        if(m_view_mode==view_FPS)//get to center
        {
            m_eye.x=m_center_pos[0];
            m_eye.y=m_center_pos[1];
            m_eye.z=m_center_pos[2];
        }
        m_eye_separation=2;
        m_convergence=30;
        m_side_shift_view=0;
        m_zoom=-50;
        init_arcball();
        m_key_delay=200;
    }
    //New settings
    if(keys[67]&&m_calc_done) //C
    {
        new_channel_edges_settings();
        keys[67]=false;
    }
    //Forced channel leakage
    if(keys[67]&&!m_void_vol_set&&m_want_void_vol) //C
    {
        m_channel_is_leaking=true;
        SetForegroundWindow(m_hWnd_console);
        keys[67]=false;
    }
    if(m_leakage_waiting_for_response)
    {
        if(keys[49]&&m_key_delay<=0) {m_leakage_option=1;m_leakage_waiting_for_response=false;} //1
        if(keys[50]&&m_key_delay<=0) {m_leakage_option=2;m_leakage_waiting_for_response=false;} //2
        if(keys[51]&&m_key_delay<=0) {m_leakage_option=3;m_leakage_waiting_for_response=false;} //3
        keys[49]=keys[50]=keys[51]=false;//reset key
    }
    //Mark block in grid
    if(m_marking_mode)
    {
        //move marker
        if(keys[100]&&m_key_delay<=0&&m_marked_block[0]<m_grid_size[0]) {m_marked_block[0]++;m_key_delay=100;}//4
        if(keys[ 97]&&m_key_delay<=0&&m_marked_block[0]>0)              {m_marked_block[0]--;m_key_delay=100;}//1
        if(keys[101]&&m_key_delay<=0&&m_marked_block[1]<m_grid_size[1]) {m_marked_block[1]++;m_key_delay=100;}//5
        if(keys[ 98]&&m_key_delay<=0&&m_marked_block[1]>0)              {m_marked_block[1]--;m_key_delay=100;}//2
        if(keys[102]&&m_key_delay<=0&&m_marked_block[2]<m_grid_size[2]) {m_marked_block[2]++;m_key_delay=100;}//6
        if(keys[ 99]&&m_key_delay<=0&&m_marked_block[2]>0)              {m_marked_block[2]--;m_key_delay=100;}//3

        if(m_channel_is_leaking) //leakage controls
        {
            //mark
            if(keys[32]&&m_key_delay<=0)//space
            {
                m_key_delay=200;
                bool marking=true;
                //remove if already marked
                for(int i=0;i<(int)m_void_blocks_marked.size();i++)
                {
                    if(m_void_blocks_marked[i].block_pos[0]==m_marked_block[0]&&
                       m_void_blocks_marked[i].block_pos[1]==m_marked_block[1]&&
                       m_void_blocks_marked[i].block_pos[2]==m_marked_block[2])
                    {
                        marking=false;
                        m_void_blocks_marked.erase(m_void_blocks_marked.begin()+i);
                        break;
                    }
                }
                if(marking)
                {
                    //add if not marked
                    st_block_pos temp_block;
                    temp_block.block_pos[0]=m_marked_block[0];
                    temp_block.block_pos[1]=m_marked_block[1];
                    temp_block.block_pos[2]=m_marked_block[2];
                    m_void_blocks_marked.push_back(temp_block);
                }
            }
            //finish
            if(keys[13]&&m_key_delay<=0)//enter
            {
                m_leakage_option=4;
                m_key_delay=200;
            }
        }
        else if(m_want_defined_channel) //channel defining
        {
            //mark
            if(keys[32]&&m_key_delay<=0&&!m_channel_center_set)//space
            {
                if(!m_channel_start_set)
                {
                    m_channel_start_pos=vec(m_marked_block[0],m_marked_block[1],m_marked_block[2]);
                    m_channel_start_set=true;
                    cout<<"Channel starting position: X = "<<m_channel_start_pos.x
                                                 <<", Y = "<<m_channel_start_pos.y
                                                 <<", Z = "<<m_channel_start_pos.z;
                    HANDLE hConsole;
                    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                    SetConsoleTextAttribute(hConsole, 14); //set YELLOW text color
                    cout<<"\nNow select the end point of the channel\n";
                    SetConsoleTextAttribute(hConsole, 15); //set WHITE text color

                }
                else if(!m_channel_end_set)
                {
                    m_channel_end_pos=vec(m_marked_block[0],m_marked_block[1],m_marked_block[2]);
                    m_channel_end_set=true;
                    cout<<"Channel ending position: X = "<<m_channel_end_pos.x
                                               <<", Y = "<<m_channel_end_pos.y
                                               <<", Z = "<<m_channel_end_pos.z;
                    HANDLE hConsole;
                    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                    SetConsoleTextAttribute(hConsole, 14); //set YELLOW text color
                    cout<<"\nNow select a position inside the channel\n";
                    SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
                }
                else if(!m_channel_center_set)
                {
                    //test position
                    if(m_pGrid[m_marked_block[0]][m_marked_block[1]][m_marked_block[2]])
                    {
                        HANDLE hConsole;
                        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
                        SetConsoleTextAttribute(hConsole, 12); //set text color
                        cout<<"\nChosen position is colliding with the protein!\n";
                        SetConsoleTextAttribute(hConsole, 15); //set text color
                    }
                    else
                    {
                        m_channel_center_pos=vec(m_marked_block[0],m_marked_block[1],m_marked_block[2]);
                        m_channel_center_set=true;
                        cout<<"A position inside the channel: X = "<<m_channel_center_pos.x
                                                         <<", Y = "<<m_channel_center_pos.y
                                                         <<", Z = "<<m_channel_center_pos.z
                            <<"\n\nPress [Enter] to finish\n";
                    }
                }
                m_key_delay=200;
            }
            //erase last marking
            if(keys[46]&&m_key_delay<=0&&m_channel_start_set)//delete
            {
                if(m_channel_center_set)
                {
                    m_channel_center_set=false;
                    cout<<"A position inside the channel cleared\n";
                }
                else if(m_channel_end_set)
                {
                    m_channel_end_set=false;
                    cout<<"Channel ending position cleared\n";
                }
                else if(m_channel_start_set)
                {
                    m_channel_start_set=false;
                    cout<<"Channel starting position cleared\n";
                }
                m_key_delay=200;
            }
            //finish
            if(keys[13]&&m_key_delay<=0&&m_channel_start_set&&m_channel_end_set&&m_channel_center_set)//enter
            {
                //test channel pos
                if(m_channel_start_pos.x>=0 && m_channel_start_pos.x<m_grid_size[0] &&
                   m_channel_start_pos.y>=0 && m_channel_start_pos.y<m_grid_size[1] &&
                   m_channel_start_pos.z>=0 && m_channel_start_pos.z<m_grid_size[2] &&
                   m_channel_end_pos.x>=0 && m_channel_end_pos.x<m_grid_size[0] &&
                   m_channel_end_pos.y>=0 && m_channel_end_pos.y<m_grid_size[1] &&
                   m_channel_end_pos.z>=0 && m_channel_end_pos.z<m_grid_size[2]) //blocks inside grid
                {
                    //test if not the same
                    if(m_channel_start_pos.x==m_channel_end_pos.x &&
                       m_channel_start_pos.y==m_channel_end_pos.y &&
                       m_channel_start_pos.z==m_channel_end_pos.z)
                    {
                        ;//should not happen
                    }
                    else//pos OK
                    {
                        //calc plane vector
                        m_channel_vector=m_channel_end_pos-m_channel_start_pos;
                        //normalize to half grid size length
                        m_channel_vector=m_channel_vector.unit()*(m_grid_gap_size/2);
                        m_channel_defined=true;
                        m_marking_mode=false;
                        cout<<"\nNew Channel Axis calculated\n";
                    }
                }
                m_key_delay=200;
            }
        }
    }

    //Scroll zoom
    switch(m_view_mode)
    {
        case view_arcball:
        {
            if(mouse_buttons[2]) m_eye.z-=2;
            if(mouse_buttons[3]) m_eye.z+=2;
        }break;
        case view_locked_roll:
        {
            if(mouse_buttons[2]) m_zoom-=2;
            if(mouse_buttons[3]) m_zoom+=2;
        }break;
    }

    //FPS movement
    if(m_view_mode==view_FPS)
    {
        //mouse drag view
        if(mouse_buttons[0])
        {
            while(m_rot_hori>360) m_rot_hori-=360;
            while(m_rot_hori<0)   m_rot_hori+=360;
            while(m_rot_vert>360) m_rot_vert-=360;
            while(m_rot_vert<0) m_rot_vert+=360;
            //rotation
            float msens=0.2;
            while(mouse_pos[0]>m_prev_mouse_pos[0])
            {
                if(m_rot_vert<180) m_rot_hori-=msens;
                else m_rot_hori+=msens;//camera is upsidedown
                m_prev_mouse_pos[0]++;
            }
            while(mouse_pos[0]<m_prev_mouse_pos[0])
            {
                if(m_rot_vert<180) m_rot_hori+=msens;
                else m_rot_hori-=msens;//camera is upsidedown
                m_prev_mouse_pos[0]--;
            }
            while(mouse_pos[1]>m_prev_mouse_pos[1])
            {
                m_rot_vert+=msens;
                m_prev_mouse_pos[1]++;
            }
            while(mouse_pos[1]<m_prev_mouse_pos[1])
            {
                m_rot_vert-=msens;
                m_prev_mouse_pos[1]--;
            }
        }
        else
        {
            m_prev_mouse_pos[0]=mouse_pos[0];
            m_prev_mouse_pos[1]=mouse_pos[1];
        }

        //key movement
        float xdir=sinf((m_rot_vert)*_piover180)*cosf((m_rot_hori)*_piover180);
        float ydir=cosf(m_rot_vert*_piover180);
        float zdir=-sinf((m_rot_vert)*_piover180)*sinf((m_rot_hori)*_piover180);
        float xrig=sinf((m_rot_vert)*_piover180)*cosf((m_rot_hori+90)*_piover180);
        float yrig=0;
        float zrig=-sinf((m_rot_vert)*_piover180)*sinf((m_rot_hori+90)*_piover180);

        float sens=0.015*cycletime;
        if(keys[87]) //w
        {
            m_eye.x+=xdir*sens;
            m_eye.y+=ydir*sens;
            m_eye.z+=zdir*sens;
        }
        if(keys[83]) //s
        {
            m_eye.x-=xdir*sens;
            m_eye.y-=ydir*sens;
            m_eye.z-=zdir*sens;
        }
        if(keys[65]) //a
        {
            m_eye.x-=xrig*sens;
            m_eye.y-=yrig*sens;
            m_eye.z-=zrig*sens;
        }
        if(keys[68]) //d
        {
            m_eye.x+=xrig*sens;
            m_eye.y+=yrig*sens;
            m_eye.z+=zrig*sens;
        }
        if(keys[81]) //q (up)
        {
            m_eye.x+=m_eye_up.x*sens;
            m_eye.y+=m_eye_up.y*sens;
            m_eye.z+=m_eye_up.z*sens;
        }
        if(keys[69]) //e (down)
        {
            m_eye.x-=m_eye_up.x*sens;
            m_eye.y-=m_eye_up.y*sens;
            m_eye.z-=m_eye_up.z*sens;
        }
    }

    //Mouse drag
    if(mouse_buttons[0] && !mouse_buttons[1]) //left MB
    {
        if(m_draging) //if draging before, move protein
        {
            switch(m_view_mode)
            {
                case view_arcball:
                {
                    m_arcball.move(m_window_length-mouse_pos[0],mouse_pos[1]);
                }break;
                case view_locked_roll:
                {
                    //Rotate right
                    while(mouse_pos[0]>m_prev_mouse_pos[0])
                    {
                        m_side_shift_view+=0.1;
                        m_prev_mouse_pos[0]+=1;
                    }
                    //Rotate left
                    while(mouse_pos[0]<m_prev_mouse_pos[0])
                    {
                        m_side_shift_view-=0.1;
                        m_prev_mouse_pos[0]-=1;
                    }
                    //Rotate up
                    while(mouse_pos[1]>m_prev_mouse_pos[1])
                    {
                        m_axis_roll+=1;
                        m_prev_mouse_pos[1]+=1;
                    }
                    //Rotate down
                    while(mouse_pos[1]<m_prev_mouse_pos[1])
                    {
                        m_axis_roll-=1;
                        m_prev_mouse_pos[1]-=1;
                    }
                }break;
                case view_FPS:
                {
                    ;
                }break;
            }
        }
        else //enable draging
        {
            m_arcball.start(m_window_length-mouse_pos[0],mouse_pos[1]);
            m_draging=true;
            //remember mouse pos
            m_prev_mouse_pos[0]=mouse_pos[0];
            m_prev_mouse_pos[1]=mouse_pos[1];
        }
    }
    else m_draging=false; //disable draging

    //Move center
    if(mouse_buttons[1] && !mouse_buttons[0])
    {
        if(m_moving)
        {
            //float sens=0.1;
            //Rotate right
            while(mouse_pos[0]>m_prev_mouse_pos[0])
            {
                m_eye.x+=0.1;
                m_prev_mouse_pos[0]+=1;
            }
            //Rotate left
            while(mouse_pos[0]<m_prev_mouse_pos[0])
            {
                m_eye.x-=0.1;
                m_prev_mouse_pos[0]-=1;
            }
            //Rotate up
            while(mouse_pos[1]>m_prev_mouse_pos[1])
            {
                m_eye.y-=0.1;
                m_prev_mouse_pos[1]+=1;
            }
            //Rotate down
            while(mouse_pos[1]<m_prev_mouse_pos[1])
            {
                m_eye.y+=0.1;
                m_prev_mouse_pos[1]-=1;
            }
        }
        else
        {
            m_moving=true;
            //remember mouse pos
            m_prev_mouse_pos[0]=mouse_pos[0];
            m_prev_mouse_pos[1]=mouse_pos[1];
        }
    }
    else m_moving=false;

    //Mouse drag zoom
    if(mouse_buttons[0] && mouse_buttons[1])
    {
        if(m_drag_zooming)
        {
            switch(m_view_mode)
            {
                case view_arcball:
                {
                    while(mouse_pos[1]<m_prev_mouse_pos[1])
                    {
                        m_eye.z-=0.2;
                        m_prev_mouse_pos[1]-=1;
                    }
                    while(mouse_pos[1]>m_prev_mouse_pos[1])
                    {
                        m_eye.z+=0.2;
                        m_prev_mouse_pos[1]+=1;
                    }
                }break;
                case view_locked_roll:
                {
                    while(mouse_pos[1]<m_prev_mouse_pos[1])
                    {
                        m_zoom-=0.2;
                        m_prev_mouse_pos[1]-=1;
                    }
                    while(mouse_pos[1]>m_prev_mouse_pos[1])
                    {
                        m_zoom+=0.2;
                        m_prev_mouse_pos[1]+=1;
                    }
                }break;
            }
        }
        else
        {
            m_drag_zooming=true;
            m_prev_mouse_pos[1]=mouse_pos[1];
        }
    }
    else m_drag_zooming=false;

    return true;
}

void protein::set_camera_pos(void)
{
    switch(m_view_mode)
    {
        case view_arcball:
        {
            //Arcball rotation
            glTranslatef(m_eye.x,m_eye.y,m_eye.z);//move to eye pos
            m_arcball.rotate();//rotate scene
            //center
            glTranslatef(-m_center_pos[0],-m_center_pos[1],-m_center_pos[2]);//move to center
            //glTranslatef(-m_focus_pos[0],-m_focus_pos[1],-m_focus_pos[2]);//move to focus
            //get camera pos
            GLfloat mdl[16];
            glGetFloatv(GL_MODELVIEW_MATRIX, mdl);
            m_cam_pos.x= -(mdl[0] * mdl[12] + mdl[1] * mdl[13] + mdl[2] * mdl[14]);
            m_cam_pos.y= -(mdl[4] * mdl[12] + mdl[5] * mdl[13] + mdl[6] * mdl[14]);
            m_cam_pos.z= -(mdl[8] * mdl[12] + mdl[9] * mdl[13] + mdl[10] * mdl[14]);
        } break;
        case view_locked_roll:
        {
            glTranslatef(m_side_shift_view,0,0); //shift view sideways
            glTranslatef(0,0,m_zoom); //zooming
            glRotatef(90,0,1,0);//side view
            if(m_want_defined_channel&&m_channel_defined)
            {
                //Normalize channel vector
                vec channel_axis=m_channel_end_pos-m_channel_start_pos;
                channel_axis=channel_axis.unit();
                //get angle between vectors (dot product between channel and z-axis)
                float rot_angle=acos(channel_axis*vec(0,0,1))/_piover180;
                //get perpendicular axis
                vec rot_axis=channel_axis^vec(0,0,1);
                rot_axis=rot_axis.unit();
                //rotate chanel axis to Z-axis
                glRotatef(rot_angle, rot_axis.x, rot_axis.y, rot_axis.z);
                //rotate around channel axis
                //view chosen angle
                glRotatef(180+m_axis_roll,channel_axis.x,channel_axis.y,channel_axis.z); //calculations done
                //get pos of center of channel
                glTranslatef(-(m_channel_start_pos.x+(m_channel_end_pos.x-m_channel_start_pos.x)*0.5)*m_grid_gap_size-m_grid_pos.x,
                             -(m_channel_start_pos.y+(m_channel_end_pos.y-m_channel_start_pos.y)*0.5)*m_grid_gap_size-m_grid_pos.y,
                             -(m_channel_start_pos.z+(m_channel_end_pos.z-m_channel_start_pos.z)*0.5)*m_grid_gap_size-m_grid_pos.z);
            }
            else
            {
                if(m_want_dist_channel&&!m_have_convex_channel_lengths&&!m_channel_is_leaking&&!m_marking_mode) glRotatef(m_axis_roll_progress,0,0,1); //chosen angle
                else                                                                                            glRotatef(m_axis_roll,0,0,1); //calculations done
                glTranslatef(-m_center_pos[0],-m_center_pos[1],-m_center_pos[2]);//move to center of protein
            }
        } break;
        case view_FPS:
        {
            //get view vectors
            m_eye_dir.x=sinf((m_rot_vert)*_piover180)*cosf((m_rot_hori)*_piover180);
            m_eye_dir.y=cosf(m_rot_vert*_piover180);
            m_eye_dir.z=-sinf((m_rot_vert)*_piover180)*sinf((m_rot_hori)*_piover180);

            m_eye_up.x=-sinf((m_rot_vert-90)*_piover180)*cosf(m_rot_hori*_piover180);
            m_eye_up.y=-cosf((m_rot_vert-90)*_piover180);
            m_eye_up.z=sinf((m_rot_vert-90)*_piover180)*sinf(m_rot_hori*_piover180);

            m_eye_right.x=sinf((m_rot_vert)*_piover180)*cosf((m_rot_hori+90)*_piover180);
            m_eye_right.y=cosf((m_rot_vert)*_piover180);
            m_eye_right.z=-sinf((m_rot_vert)*_piover180)*sinf((m_rot_hori+90)*_piover180);

            gluLookAt(//move camera
                      m_eye.x,m_eye.y,m_eye.z,
                      m_eye.x+m_eye_dir.x,m_eye.y+m_eye_dir.y,m_eye.z+m_eye_dir.z,
                      m_eye_up.x,m_eye_up.y,m_eye_up.z);
        } break;
    }
}

atom protein::get_atom_from_string(string line)
{
    //Send string line data to atom object
    atom temp_atom;
    //Get coord
    string sNumber;
    sNumber.append(1,line[31]);
    sNumber.append(1,line[32]);
    sNumber.append(1,line[33]);
    sNumber.append(1,line[34]);
    sNumber.append(1,line[35]);
    sNumber.append(1,line[36]);
    sNumber.append(1,line[37]);
    temp_atom.x=atof(sNumber.c_str());
    sNumber.clear();
    sNumber.append(1,line[39]);
    sNumber.append(1,line[40]);
    sNumber.append(1,line[41]);
    sNumber.append(1,line[42]);
    sNumber.append(1,line[43]);
    sNumber.append(1,line[44]);
    sNumber.append(1,line[45]);
    temp_atom.y=atof(sNumber.c_str());
    sNumber.clear();
    sNumber.append(1,line[47]);
    sNumber.append(1,line[48]);
    sNumber.append(1,line[49]);
    sNumber.append(1,line[50]);
    sNumber.append(1,line[51]);
    sNumber.append(1,line[52]);
    sNumber.append(1,line[53]);
    temp_atom.z=atof(sNumber.c_str());
    //Get radius
    temp_atom.type=line[13];
    float radius;
    switch (temp_atom.type)
    {
        case 'C': radius=0.7; break;
        case 'O': radius=0.6; break;
        case 'N': radius=0.65;break;
        case 'S': radius=1.0; break;
        case 'H': radius=0.25;break;
        default: radius=1;
    }
    temp_atom.radius=radius;
    return temp_atom;
}

bool protein::collision_water(float xpos,float ypos,float zpos)
{
    //check if a water molecule could be placed at given point without colliding with the protein
    float solvent_radius=m_solvent_radius;
    for(int atom=0;atom<(int)m_atoms.size();atom++)
    {
        //quick test (treat as cubes)
        if(m_atoms[atom].x>xpos+m_atoms[atom].radius+solvent_radius || m_atoms[atom].x<xpos-m_atoms[atom].radius-solvent_radius ||
           m_atoms[atom].y>ypos+m_atoms[atom].radius+solvent_radius || m_atoms[atom].y<ypos-m_atoms[atom].radius-solvent_radius ||
           m_atoms[atom].z>zpos+m_atoms[atom].radius+solvent_radius || m_atoms[atom].z<zpos-m_atoms[atom].radius-solvent_radius)
        {
            continue;
        }
        else//treat as spheres
        {
            vec temp_vec(m_atoms[atom].x-xpos,m_atoms[atom].y-ypos,m_atoms[atom].z-zpos);
            if(temp_vec.length()<m_atoms[atom].radius+solvent_radius) return true;
        }
    }
    return false; //no collision
}

void protein::init_arcball(void)
{
    m_arcball.reset();
    float zoom=0.05;//rotation sensitivity
    m_eye=vec(0,0,-50);
    m_center=vec(0,0,0);
    m_up=vec(0,1,0);
    m_arcball.setzoom(zoom,m_eye,m_up);
}

bool protein::loadMesh(const char* file_name)
{
    //call CGAL
    string mesh_gen_call("data\\MeshGenerator.exe ");
    mesh_gen_call.append(file_name);//give name of protein file
    if(system(mesh_gen_call.c_str())!=0)//get CGAL data
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR loading mesh, file not found\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //Read CGAL data
    string line;
    ifstream file("mesh_output.txt",ios::binary);
    if (!file)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR loading mesh, file not found\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //read header
    getline(file,line);//read OFF
    getline(file,line);//read number of vertices and triangles
    istringstream test_iss(line);
    string temp_value;
    test_iss>>temp_value;
    int number_of_vertices=atoi(temp_value.c_str());
    test_iss>>temp_value;
    int number_of_triangles=atoi(temp_value.c_str());
    getline(file,line);//skip empty line
    //load vertices
    for(int counter=0;counter<number_of_vertices;counter++)
    {
        getline(file,line);
        vec temp_vec;
        istringstream iss(line);
        string sFloat;
        iss>>sFloat;
        temp_vec.x=atof(sFloat.c_str());
        iss>>sFloat;
        temp_vec.y=atof(sFloat.c_str());
        iss>>sFloat;
        temp_vec.z=atof(sFloat.c_str());
        m_mesh_atom_vertices.push_back(temp_vec);//store vertex
    }
    cout<<"Protein mesh data:\n";
    cout<<"Number of vertices: "<<m_mesh_atom_vertices.size()<<endl;
    //load triangles
    for(int counter=0;counter<number_of_triangles;counter++)
    {
        getline(file,line);
        st_triangle temp_triangle;
        istringstream iss(line);
        string sIndex;
        iss>>sIndex;//skip 3
        iss>>sIndex;
        temp_triangle.first=atoi(sIndex.c_str());
        iss>>sIndex;
        temp_triangle.second=atoi(sIndex.c_str());
        iss>>sIndex;
        temp_triangle.third=atoi(sIndex.c_str());
        m_mesh_atom_triangles.push_back(temp_triangle);//store triangle
    }
    cout<<"Number of triangles: "<<m_mesh_atom_triangles.size()<<endl;
    file.close();

    return true;
}

bool protein::load_mesh_channel(const char* file_name)
{
    //call CGAL
    string mesh_gen_call("data\\MeshGenerator.exe data\\");
    mesh_gen_call.append(file_name);//give name of protein file
    if(system(mesh_gen_call.c_str())!=0)//get CGAL data
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR loading channel mesh, file not found\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //Read CGAL data
    string line;
    ifstream file("mesh_output.txt",ios::binary);
    if (!file)
    {
        HANDLE hConsole;
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
        SetConsoleTextAttribute(hConsole, 12); //set RED text color
        cout<<"ERROR loading channel mesh, could not create file\n";
        SetConsoleTextAttribute(hConsole, 15); //set WHITE text color
        return false;
    }
    //read header
    getline(file,line);//read OFF
    getline(file,line);//read number of vertices and triangles
    istringstream test_iss(line);
    string temp_value;
    test_iss>>temp_value;
    int number_of_vertices=atoi(temp_value.c_str());
    test_iss>>temp_value;
    int number_of_triangles=atoi(temp_value.c_str());
    getline(file,line);//skip empty line
    //load vertices
    for(int counter=0;counter<number_of_vertices;counter++)
    {
        getline(file,line);
        vec temp_vec;
        istringstream iss(line);
        string sFloat;
        iss>>sFloat;
        temp_vec.x=atof(sFloat.c_str());
        iss>>sFloat;
        temp_vec.y=atof(sFloat.c_str());
        iss>>sFloat;
        temp_vec.z=atof(sFloat.c_str());
        m_mesh_channel_vertices.push_back(temp_vec);//store vertex
    }
    cout<<"Channel mesh data:\n";
    cout<<"Number of vertices: "<<m_mesh_channel_vertices.size()<<endl;
    //load triangles
    for(int counter=0;counter<number_of_triangles;counter++)
    {
        getline(file,line);
        st_triangle temp_triangle;
        istringstream iss(line);
        string sIndex;
        iss>>sIndex;//skip 3
        iss>>sIndex;
        temp_triangle.first=atoi(sIndex.c_str());
        iss>>sIndex;
        temp_triangle.second=atoi(sIndex.c_str());
        iss>>sIndex;
        temp_triangle.third=atoi(sIndex.c_str());
        m_mesh_channel_triangles.push_back(temp_triangle);//store triangle
    }
    cout<<"Number of triangles: "<<m_mesh_channel_triangles.size()<<endl;
    file.close();

    return true;
}

bool protein::new_channel_edges_settings(void)
{
    SetForegroundWindow(m_hWnd_console); //make console window active
    m_calc_done=false;
    //reset grid
    m_void_blocks.clear();
    m_void_blocks_exp.clear();
    m_void_blocks_new.clear();
    m_void_blocks_marked.clear();
    delete[] m_pGrid;
    delete[] m_pGrid_color;
    m_pGrid=0;
    m_pGrid_color=0;
    m_grid_updated=false;
    m_have_grid_void_color=false;
    m_grid_z_progress=0;
    m_protein_volume=0;
    //reset void data
    m_void_vol_set=false;
    m_void_vol_progress=0;
    //reset channel edges data
    m_tried_to_get_channel_mesh=false;
    m_channel_edge_set=false;
    m_calc_done=false;
    m_channel_edges_up.clear();
    m_channel_edges_down.clear();
    m_channel_edges_normals.clear();
    m_channel_edges_color.clear();
    m_channel_edges_length.clear();
    m_axis_roll_progress=0;
    m_have_convex_channel_lengths=false;
    //reset channel mesh
    m_mesh_channel_vertices.clear();
    m_mesh_channel_triangles.clear();
    m_vertex_triangle_map_channel.clear();
    m_have_mesh_channel_normals=false;
    m_have_real_mesh_channel_normals=false;
    m_mesh_normal_progress=0;
    m_mesh_normal_vertex_progress=0;
    //reset channel placement
    m_channel_defined=false;
    m_channel_start_set=false;
    m_channel_end_set=false;
    m_channel_center_set=false;

    //get new settings
    float grid_gap_size=0;
    float exp_size=0;
    float step_size=0;
    float solvent_radius=0;
    int   exposure_lvl=0;
    string input;

    HANDLE hConsole;
    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
    SetConsoleTextAttribute(hConsole, 14); //set text color
    cout<<"\nEnter value for grid box size in Angstom (default 1.0) : ";
    SetConsoleTextAttribute(hConsole, 15); //set text color
    getline(cin,input);
    SetConsoleTextAttribute(hConsole, 14); //set text color
    grid_gap_size=atof(input.c_str());
    if(m_want_void_vol&&m_want_defined_channel)
    {
        cout<<"\nDescription of exposure level cut-off for channel/cleft detection:\n"
              "The exposure level is how many walls of the grid box that the current\n"
              "grid position can \"see\" without passing through parts of the protein.\n"
              "Grid positions with higher exposure level then selected will be ignored.\n\n";
        cout<<"Enter value for exposeure level cut-off for user defined channel (0-6) : ";
        SetConsoleTextAttribute(hConsole, 15); //set text color
        getline(cin,input);
        SetConsoleTextAttribute(hConsole, 14); //set text color
        exposure_lvl=atoi(input.c_str());
    }
    cout<<"Enter value for solvent radius (default 1.4 for water) : ";
    SetConsoleTextAttribute(hConsole, 15); //set text color
    getline(cin,input);
    SetConsoleTextAttribute(hConsole, 14); //set text color
    solvent_radius=atof(input.c_str());
    if(m_want_dist_channel)
    {
        cout<<"Enter channel ruler step size along channel axis (default 0.5) : ";
        SetConsoleTextAttribute(hConsole, 15); //set text color
        getline(cin,input);
        SetConsoleTextAttribute(hConsole, 14); //set text color
        step_size=atof(input.c_str());
        cout<<"Enter channel ruler side expansion size (default 0.5) : ";
        SetConsoleTextAttribute(hConsole, 15); //set text color
        getline(cin,input);
        exp_size=atof(input.c_str());
    }
    cout<<"\nRecalculating...\n\n";

    //new values must be positive
    if(grid_gap_size>0.01) m_grid_gap_size=grid_gap_size;
    if(exp_size>0) m_channel_edges_exp_size=exp_size;
    if(step_size>0) m_channel_edges_step_size=step_size;
    if(solvent_radius>0) m_solvent_radius=solvent_radius;
    if(m_want_defined_channel)
    {
        if(exposure_lvl<7&&exposure_lvl>=0) m_show_grid_void_exposure_level=exposure_lvl;
    }

    SetForegroundWindow(m_hWnd_gui);
    return true;
}

void protein::setWindowVal(int winWidth,int winHeight,int scrWidth,int scrHeight,HWND hwnd_console,HWND hwnd_gui)
{
    m_hWnd_console=hwnd_console;
    m_hWnd_gui=hwnd_gui;
    m_window_length=winWidth;
    m_window_height=winHeight;
    m_screen_length=scrWidth;
    m_screen_height=scrHeight;
    init_arcball(); //Arcball needs to be reseted
    Text_box.reload(); //textures needs to be reloaded
}

void protein::setWindowVal(int winWidth,int winHeight,int scrWidth,int scrHeight)
{
    m_window_length=winWidth;
    m_window_height=winHeight;
    m_screen_length=scrWidth;
    m_screen_height=scrHeight;
    init_arcball(); //Arcball needs to be reseted
    Text_box.reload(); //textures needs to be reloaded
}

bool protein::is_block_inside_planes(st_block_pos block)
{
    float distA,distB;//relative distances without sqrt(), only for comparison
    //check first plane
    distA=(m_channel_start_pos.x-(float)block.block_pos[0])*(m_channel_start_pos.x-(float)block.block_pos[0])+
          (m_channel_start_pos.y-(float)block.block_pos[1])*(m_channel_start_pos.y-(float)block.block_pos[1])+
          (m_channel_start_pos.z-(float)block.block_pos[2])*(m_channel_start_pos.z-(float)block.block_pos[2]);
    distB=(m_channel_start_pos.x+m_channel_vector.x-(float)block.block_pos[0])*(m_channel_start_pos.x+m_channel_vector.x-(float)block.block_pos[0])+
          (m_channel_start_pos.y+m_channel_vector.y-(float)block.block_pos[1])*(m_channel_start_pos.y+m_channel_vector.y-(float)block.block_pos[1])+
          (m_channel_start_pos.z+m_channel_vector.z-(float)block.block_pos[2])*(m_channel_start_pos.z+m_channel_vector.z-(float)block.block_pos[2]);
    if(distA<distB) return false;//wrong side of plane
    //check other plane
    distA=(m_channel_end_pos.x-(float)block.block_pos[0])*(m_channel_end_pos.x-(float)block.block_pos[0])+
          (m_channel_end_pos.y-(float)block.block_pos[1])*(m_channel_end_pos.y-(float)block.block_pos[1])+
          (m_channel_end_pos.z-(float)block.block_pos[2])*(m_channel_end_pos.z-(float)block.block_pos[2]);
    distB=(m_channel_end_pos.x-m_channel_vector.x-(float)block.block_pos[0])*(m_channel_end_pos.x-m_channel_vector.x-(float)block.block_pos[0])+
          (m_channel_end_pos.y-m_channel_vector.y-(float)block.block_pos[1])*(m_channel_end_pos.y-m_channel_vector.y-(float)block.block_pos[1])+
          (m_channel_end_pos.z-m_channel_vector.z-(float)block.block_pos[2])*(m_channel_end_pos.z-m_channel_vector.z-(float)block.block_pos[2]);
    if(distA<distB) return false;//wrong side of plane

    return true;//block inside planes
}



vec rotatePointAboutLine(vec p,float theta,vec p1,vec p2)
{
   //Rotation of a point in 3 dimensional space by theta about an arbitrary axes defined by a line
   //between two points P1 = (x1,y1,z1) and P2 = (x2,y2,z2) can be achieved by the following steps
   vec u,q1,q2;
   float d;

   // Step 1 translate space so that the rotation axis passes through the origin
   q1.x = p.x - p1.x;
   q1.y = p.y - p1.y;
   q1.z = p.z - p1.z;

   u.x = p2.x - p1.x;
   u.y = p2.y - p1.y;
   u.z = p2.z - p1.z;
   //Normalise(&u);
   u=u.unit();
   d = sqrt(u.y*u.y + u.z*u.z);

   // Step 2 rotate space about the x axis so that the rotation axis lies in the xz plane
   if (d != 0)
   {
      q2.x = q1.x;
      q2.y = q1.y * u.z / d - q1.z * u.y / d;
      q2.z = q1.y * u.y / d + q1.z * u.z / d;
   }
   else
   {
      q2 = q1;
   }

   // Step 3 rotate space about the y axis so that the rotation axis lies along the z axis
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;

   // Step 4 perform the desired rotation by theta about the z axis
   q2.x = q1.x * cosf(theta) - q1.y * sinf(theta);
   q2.y = q1.x * sinf(theta) + q1.y * cosf(theta);
   q2.z = q1.z;

   // Inverse of step 3
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   // Inverse of step 2
   if (d != 0)
   {
      q2.x =   q1.x;
      q2.y =   q1.y * u.z / d + q1.z * u.y / d;
      q2.z = - q1.y * u.y / d + q1.z * u.z / d;
   }
   else
   {
      q2 = q1;
   }

   // Inverse of step 1
   q1.x = q2.x + p1.x;
   q1.y = q2.y + p1.y;
   q1.z = q2.z + p1.z;

   return q1;
}

float volume_of_triangle(vec p1, vec p2, vec p3)
{
    float v321 = p3.x*p2.y*p1.z;
    float v231 = p2.x*p3.y*p1.z;
    float v312 = p3.x*p1.y*p2.z;
    float v132 = p1.x*p3.y*p2.z;
    float v213 = p2.x*p1.y*p3.z;
    float v123 = p1.x*p2.y*p3.z;
    return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

bool float_to_char(float value,char cDigits[8])
{
    value+=0.0005;//for rounding
    int iDigits[7];
    if(value>99||value<-99) return false;//too big value;
    bool is_neg=false;
    if(value<0) is_neg=true;
    float multi=0.1;
    for(int i=0;i<7;i++)
    {
        iDigits[i]=abs(int(value*multi)-int(value*multi/10)*10);
        multi*=10;
    }
    //convert to char
    cDigits[1]=int_2_char(iDigits[0]);
    cDigits[2]=int_2_char(iDigits[1]);
    cDigits[4]=int_2_char(iDigits[2]);
    cDigits[5]=int_2_char(iDigits[3]);
    cDigits[6]=int_2_char(iDigits[4]);
    if(is_neg)
    {
        //place minus sign
        if(iDigits[0]>0) cDigits[0]='-';
        else
        {
            cDigits[0]=' ';
            cDigits[1]='-';
        }
    }
    else
    {
        if(iDigits[0]>0) cDigits[0]=' ';
        else
        {
            cDigits[0]=' ';
            cDigits[1]=' ';
        }
    }
    return true;
}

char int_2_char(int value)
{
    switch (value)
    {
        case 0: return '0';
        case 1: return '1';
        case 2: return '2';
        case 3: return '3';
        case 4: return '4';
        case 5: return '5';
        case 6: return '6';
        case 7: return '7';
        case 8: return '8';
        case 9: return '9';
    }
    return 'X';//error
}

void get_color(int value,float color[3])
{
    //Get a random color that will be based on the int value given
    srand(value);
    color[0]=float(rand()%100)/100;
    color[1]=float(rand()%100)/100;
    color[2]=float(rand()%100)/100;
}

vec get_3D_pos(int x, int y)
{
    //get point where coordinates hit an object in the window (not used right now)
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLfloat winX, winY, winZ;
    GLdouble posX, posY, posZ;

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );

    winX = (float)x;
    winY = (float)viewport[3] - (float)y;
    glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );

    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
    return vec(posX, posY, posZ);
}

void print_controls(void)
{
    HANDLE hConsole;
    hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
    SetConsoleTextAttribute(hConsole, 10); //set text color
    cout<<"\n"
    "Loading is Complete\n";
    SetConsoleTextAttribute(hConsole, 15); //set text color
    cout<<""
    "*******************\n"
    "Controls:\n"
    "[F1]        - Fullscreen\n"
    "[F2]        - Anaglyph 3D\n"
    "[F3]        - Interlaced 3D\n"
    "[F4]        - Swap Eye Position For 3D\n"
    "[Arrows]    - Adjust Eye Separation (Right/Left) and Convergence (Up/Down)\n"
    "[Mouse]     - Rotate Protein\n"
    "[Scroll]    - Zooming\n"
    "[Backspace] - Reset View\n"
    "[Tab]       - Change View Mode\n"
    "[Num+/-]    - Rotate Channel Ruler\n"
    "[PgUp/Dn]   - Increase/Decrease Level of Exposure Tolerance\n"
    "[R]         - Load a new Protein\n"
    "[Esc]       - Exit Program\n"
    "[B]         - Toggle Between Black and White Background\n"
    "[T]         - Toggle Mesh Transparency\n"
    "[G]         - Toogle See-through Channel Mesh\n"
    "[1]         - Show Atoms\n"
    "[2]         - Show Alpha Carbon Chains\n"
    "[3]         - Show Protein Grid\n"
    "[4]         - Show Channel Grid\n"
    "[5]         - Show Channel Radiuses\n"
    "[6]         - Show Channel Ruler\n"
    "[7]         - Show Protein Mesh\n"
    "[8]         - Show Channel Mesh\n"
    "[9]         - Show Level of Exposure for Clefts\n"
    "[0]         - Show Other Atoms\n"
    "[L]         - Show Amino Acid Labels\n"
    "[F]         - Define Channel Position\n"
    "[C]         - Change Settings\n"
    "[P]         - Output Diameters of the Channel for the Current Ruler Axis\n"
    "[O]         - Output Average Diameters of Channel Along the Channel\n"
    "[I]         - Output Minimum Diameters of Channel Along the Channel\n"
    "[U]         - Output Minimum Distances from the Center Axis to the Channel Edge\n"
    "*******************\n\n";
}
