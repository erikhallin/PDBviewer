//This software allows loading of protein PDB files,
//viewing and measuring size and shape of protein channels.
//Created by Erik Hallin 2012

#define _WIN32_WINNT 0x0500 // Needed for GetConsoleWindow()
#include <windows.h>		// Header File For Windows
#include <fcntl.h>          // Needed for AllocConsole()
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include "protein.h"
#include "menu.h"
#include "resource.h"

using namespace std;

HDC			hDC=NULL;		   // Private GDI Device Context
HGLRC		hRC=NULL;		   // Permanent Rendering Context
HWND		hWnd=NULL;		   // Holds Our GUI Window Handle
HWND        hWnd_console=NULL; // Holds Our Console Handle
HINSTANCE	hInstance;		   // Holds The Instance Of The Application

//Input variables
static bool	g_keys[256];
static int  g_mouse_pos[2];
static bool g_mouse_buttons[4];
int         g_clickDelay=0;
//Window
bool	active=TRUE;		// Window Active Flag Set To TRUE By Default
bool	fullscreen=FALSE;	// Fullscreen Flag Set To Fullscreen Mode By Default
int     windowWidth=640;
int     windowHeight=480;
int     screenWidth=640;
int     screenHeight=480;
BOOL	g_exit=FALSE;	    // Bool Variable To Exit Loop
//Timing
DWORD iTickCount=0;
DWORD iLastTickCount;
float FPS=0;
float g_cycleTime=1;
//==================================
enum     states{stSTARTUP,stINPUT,stLOADING,stRUNNING,stEXIT};
int      g_state=stSTARTUP;
protein* pProtein;
menu     menu_get_file_name;
//==================================

//Function Def
void draw_input_menu(void);
void keyCheck(void);
LRESULT	CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// Declaration For WndProc

GLvoid ReSizeGLScene(GLsizei width, GLsizei height)		// Resize And Initialize The GL Window
{
	if (height==0)										// Prevent A Divide By Zero By
	{
		height=1;										// Making Height Equal One
	}
	glViewport(0,0,width,height);						// Reset The Current Viewport
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix
	gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,200.0f);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
	windowWidth=width;
	windowHeight=height;
	if(pProtein)pProtein->setWindowVal(windowWidth,windowHeight,screenWidth,screenHeight);// Tell Protein about window size
}

int InitGL(void)										// All Setup For OpenGL Goes Here
{
    screenWidth=GetSystemMetrics(SM_CXSCREEN); // Get Width of Screen
    screenHeight=GetSystemMetrics(SM_CYSCREEN);// Get Height of Screen
    //FOG
    GLfloat fogColor[4]= {0.01f, 0.01f, 0.01f, 1.0f};   // Fog Color
    glFogi(GL_FOG_MODE, GL_LINEAR);                     // Fog Mode GL_LINEAR
    glFogfv(GL_FOG_COLOR, fogColor);                    // Set Fog Color
    glFogf(GL_FOG_DENSITY, 0);                          // How Dense Will The Fog Be
    glFogf(GL_FOG_START, 10.0f);                        // Fog Start Depth
    glFogf(GL_FOG_END, 100.0f);                         // Fog End Depth
    glEnable(GL_FOG);                                   // Enables GL_FOG
    //Blending
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //Lighting
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    GLfloat LightAmbient[]=		{ 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat LightDiffuse[]=		{ 0.7f, 0.7f, 0.7f, 1.0f };
    GLfloat LightSpecular[]=    { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat LightPosition[]=	{ 0.0f, 0.0f, 0.0f, 1.0f };
	glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient);		// Setup The Ambient Light
	glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);		// Setup The Diffuse Light
	glLightfv(GL_LIGHT0, GL_SPECULAR,LightSpecular);	// Setup The Diffuse Light
	glLightfv(GL_LIGHT0, GL_POSITION,LightPosition);	// Position The Light
	glEnable(GL_LIGHT0);								// Enable Light One

	return TRUE;
}

int ProgramCycle(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer
    glLoadIdentity();

    switch (g_state)
    {
        case stSTARTUP:
        {
            cout<<//Starting message
            "\n**********************\n"
            "Protein Channel Viewer\n"
            "**********************\n"
            "This software allow viewing of protein channels, calculate the channel's volume and measuring distances inside the channel.\n"
            "The input file must be a pdb file. If the protein channel position is not\ndefined it will be assumed to be located along the Z-axis.\n\n"
            "Created by Erik Hallin\n\n";
            g_state=stINPUT;
        }break;
        case stINPUT:
        {
            HANDLE hConsole;
            hConsole = GetStdHandle(STD_OUTPUT_HANDLE); //get console handle;
            SetConsoleTextAttribute(hConsole, 14); //set text color
            cout<<"Enter name of input file located inside the \"INPUT\" folder\n(test.pdb for testing)\n\n"
                  "If not all measurments are required, specify the wanted features by writing\n\"1\" for wanted or \"0\" for unwanted, after the name of the file.\n"
                  "\nLoading Options:\n"
                  "Mesh for Atoms      - Surface for protein could be displayed\n"
                  "Mesh for Channel    - Surface of protein channel could be displayed\n"
                  "Volume of Channel   - Volume of the channel will be calculated\n"
                  "Radiuses of Channel - Distance inside the channel will be calculated\n"
                  "Channel/Cleft Axis  - User define location and axis of channel\n\n"
                  "Example: \"test.pdb 0 1 1 0 0\" will load the channel mesh and calculate the\n"
                  "volume of the channel but will not give the atom mesh nor the radiuses of the\n"
                  "channel. And the channel will use the Z-axis as channel axis.\n\n"
                  "Input file: ";
            SetConsoleTextAttribute(hConsole, 15); //set text color
            string path;
            getline(cin,path);
            string folder("INPUT\\");
            folder.append(path);
            if (menu_get_file_name.set_filename(folder)) g_state=stLOADING; //File exists
        }break;
        case stLOADING:
        {
            string filename=menu_get_file_name.get_filename();
            pProtein=new protein(filename,
                                 hWnd,
                                 hWnd_console,
                                 menu_get_file_name.m_want_mesh_atom,
                                 menu_get_file_name.m_want_mesh_channel,
                                 menu_get_file_name.m_want_void_vol,
                                 menu_get_file_name.m_want_dist_channel,
                                 menu_get_file_name.m_want_defined_channel
                                 );//create protein object
            pProtein->setWindowVal(windowWidth,windowHeight,screenWidth,screenHeight);
            if (pProtein)
            {
                g_state=stRUNNING;
            }
            else g_state=stEXIT;
        }break;
        case stRUNNING:
        {
            pProtein->update_protein();
        }break;
        case stEXIT:
        {
            delete pProtein;
            g_exit=true;
        }break;
    }

	return TRUE;
}

//GLOBAL FUNCTIONS

void keyCheck(void)
{
    switch(g_state)
    {
        case stINPUT:
        {
            //Nothing
        }break;
        case stRUNNING:
        {
            if(pProtein)
            {
                if(!pProtein->move(g_keys,g_mouse_buttons,g_mouse_pos,g_cycleTime)) //Move protein
                {
                    //load a new protein
                    //pProtein->~protein(); //destroy old protein
                    delete pProtein;
                    g_state=stINPUT;
                    g_keys[82]=false; //reset reload key 'R'
                }
            }
        }break;
    }

    g_mouse_buttons[2]=false;g_mouse_buttons[3]=false; //reset scroll
}

//==========================================================================================================================
GLvoid KillGLWindow(void)								// Properly Kill The Window
{
	if (fullscreen)										// Are We In Fullscreen Mode?
	{
		ChangeDisplaySettings(NULL,0);					// If So Switch Back To The Desktop
		ShowCursor(TRUE);								// Show Mouse Pointer
	}

	if (hRC)											// Do We Have A Rendering Context?
	{
		if (!wglMakeCurrent(NULL,NULL))					// Are We Able To Release The DC And RC Contexts?
		{
			MessageBox(NULL,"Release Of DC And RC Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}

		if (!wglDeleteContext(hRC))						// Are We Able To Delete The RC?
		{
			MessageBox(NULL,"Release Rendering Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}
		hRC=NULL;										// Set RC To NULL
	}

	if (hDC && !ReleaseDC(hWnd,hDC))					// Are We Able To Release The DC
	{
		MessageBox(NULL,"Release Device Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hDC=NULL;										// Set DC To NULL
	}

	if (hWnd && !DestroyWindow(hWnd))					// Are We Able To Destroy The Window?
	{
		MessageBox(NULL,"Could Not Release hWnd.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hWnd=NULL;										// Set hWnd To NULL
	}

	if (!UnregisterClass("Protein Channel Viewer",hInstance))			// Are We Able To Unregister Class
	{
		MessageBox(NULL,"Could Not Unregister Class.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hInstance=NULL;									// Set hInstance To NULL
	}
}

/*	This Code Creates Our OpenGL Window.  Parameters Are:					*
 *	title			- Title To Appear At The Top Of The Window				*
 *	width			- Width Of The GL Window Or Fullscreen Mode				*
 *	height			- Height Of The GL Window Or Fullscreen Mode			*
 *	bits			- Number Of Bits To Use For Color (8/16/24/32)			*
 *	fullscreenflag	- Use Fullscreen Mode (TRUE) Or Windowed Mode (FALSE)	*/

BOOL CreateGLWindow(char* title, int width, int height, int bits, bool fullscreenflag)
{
	GLuint		PixelFormat;			// Holds The Results After Searching For A Match
	WNDCLASS	wc;						// Windows Class Structure
	DWORD		dwExStyle;				// Window Extended Style
	DWORD		dwStyle;				// Window Style
	RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left=(long)0;			// Set Left Value To 0
	WindowRect.right=(long)width;		// Set Right Value To Requested Width
	WindowRect.top=(long)0;				// Set Top Value To 0
	WindowRect.bottom=(long)height;		// Set Bottom Value To Requested Height
	fullscreen=fullscreenflag;			// Set The Global Fullscreen Flag

	hInstance			= GetModuleHandle(NULL);				// Grab An Instance For Our Window
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc		= (WNDPROC) WndProc;					// WndProc Handles Messages
	wc.cbClsExtra		= 0;									// No Extra Window Data
	wc.cbWndExtra		= 0;									// No Extra Window Data
	wc.hInstance		= hInstance;							// Set The Instance
	wc.hIcon            = LoadIcon(hInstance,MAKEINTRESOURCE(IDI_PCV));// Load The Icon
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
	wc.hbrBackground	= NULL;									// No Background Required For GL
	wc.lpszMenuName		= NULL;									// We Don't Want A Menu
	wc.lpszClassName	= "Protein Channel Viewer";								// Set The Class Name

	if (!RegisterClass(&wc))									// Attempt To Register The Window Class
	{
		MessageBox(NULL,"Failed To Register The Window Class.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;											// Return FALSE
	}

	if (fullscreen)												// Attempt Fullscreen Mode?
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings,0,sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
		dmScreenSettings.dmSize=sizeof(dmScreenSettings);		// Size Of The Devmode Structure
		dmScreenSettings.dmPelsWidth	= width;				// Selected Screen Width
		dmScreenSettings.dmPelsHeight	= height;				// Selected Screen Height
		dmScreenSettings.dmBitsPerPel	= bits;					// Selected Bits Per Pixel
		dmScreenSettings.dmFields=DM_BITSPERPEL|DM_PELSWIDTH|DM_PELSHEIGHT;

		// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
		if (ChangeDisplaySettings(&dmScreenSettings,CDS_FULLSCREEN)!=DISP_CHANGE_SUCCESSFUL)
		{
			// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
			if (MessageBox(NULL,"The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?","NeHe GL",MB_YESNO|MB_ICONEXCLAMATION)==IDYES)
			{
				fullscreen=FALSE;		// Windowed Mode Selected.  Fullscreen = FALSE
			}
			else
			{
				// Pop Up A Message Box Letting User Know The Program Is Closing.
				MessageBox(NULL,"Program Will Now Close.","ERROR",MB_OK|MB_ICONSTOP);
				return FALSE;									// Return FALSE
			}
		}
	}

	if (fullscreen)												// Are We Still In Fullscreen Mode?
	{
		dwExStyle=WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle=WS_POPUP;										// Windows Style
	}
	else
	{
		dwExStyle=WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
		dwStyle=WS_OVERLAPPEDWINDOW;							// Windows Style
	}

	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);		// Adjust Window To True Requested Size

	// Create The Window
	if (!(hWnd=CreateWindowEx(	dwExStyle,							// Extended Style For The Window
								"Protein Channel Viewer",							// Class Name
								title,								// Window Title
								dwStyle |							// Defined Window Style
								WS_CLIPSIBLINGS |					// Required Window Style
								WS_CLIPCHILDREN,					// Required Window Style
								0, 0,								// Window Position
								WindowRect.right-WindowRect.left,	// Calculate Window Width
								WindowRect.bottom-WindowRect.top,	// Calculate Window Height
								NULL,								// No Parent Window
								NULL,								// No Menu
								hInstance,							// Instance
								NULL)))								// Dont Pass Anything To WM_CREATE
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Window Creation Error.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	static	PIXELFORMATDESCRIPTOR pfd=				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		bits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	if (!(hDC=GetDC(hWnd)))							// Did We Get A Device Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Device Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(PixelFormat=ChoosePixelFormat(hDC,&pfd)))	// Did Windows Find A Matching Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Find A Suitable PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!SetPixelFormat(hDC,PixelFormat,&pfd))		// Are We Able To Set The Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Set The PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(hRC=wglCreateContext(hDC)))				// Are We Able To Get A Rendering Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!wglMakeCurrent(hDC,hRC))					// Try To Activate The Rendering Context
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Activate The GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	ShowWindow(hWnd,SW_SHOW);						// Show The Window
	SetForegroundWindow(hWnd);						// Slightly Higher Priority
	SetFocus(hWnd);									// Sets Keyboard Focus To The Window
	ReSizeGLScene(width, height);					// Set Up Our Perspective GL Screen

	if (!InitGL())									// Initialize Our Newly Created GL Window
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Initialization Failed.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	return TRUE;									// Success
}

LRESULT CALLBACK WndProc(	HWND	hWnd,			// Handle For This Window
							UINT	uMsg,			// Message For This Window
							WPARAM	wParam,			// Additional Message Information
							LPARAM	lParam)			// Additional Message Information
{
	switch (uMsg)									// Check For Windows Messages
	{
		case WM_ACTIVATE:							// Watch For Window Activate Message
		{
			if (!HIWORD(wParam))					// Check Minimization State
			{
				active=TRUE;						// Program Is Active
			}
			else
			{
				active=FALSE;						// Program Is No Longer Active
			}

			return 0;								// Return To The Message Loop
		}

		case WM_SYSCOMMAND:							// Intercept System Commands
		{
			switch (wParam)							// Check System Calls
			{
				case SC_SCREENSAVE:					// Screensaver Trying To Start?
				case SC_MONITORPOWER:				// Monitor Trying To Enter Powersave?
				return 0;							// Prevent From Happening
			}
			break;									// Exit
		}

		case WM_CLOSE:								// Did We Receive A Close Message?
		{
			PostQuitMessage(0);						// Send A Quit Message
			return 0;								// Jump Back
		}

        case WM_MOUSEMOVE:
        {
             g_mouse_pos[0]=LOWORD(lParam);
             g_mouse_pos[1]=HIWORD(lParam);
             return 0;
        }

        case WM_LBUTTONDOWN:
        {
             g_mouse_buttons[0]=true;
             return 0;
        }

        case WM_RBUTTONDOWN:
        {
             g_mouse_buttons[1]=true;
             return 0;
        }

        case WM_LBUTTONUP:
        {
             g_mouse_buttons[0]=false;
             return 0;
        }

        case WM_RBUTTONUP:
        {
             g_mouse_buttons[1]=false;
             return 0;
        }

        case WM_MOUSEWHEEL:
        {
            if(HIWORD(wParam)>5000) {g_mouse_buttons[2]=true;}
            if(HIWORD(wParam)>100&&HIWORD(wParam)<5000) {g_mouse_buttons[3]=true;}
            return 0;
        }

		case WM_KEYDOWN:							// Is A Key Being Held Down?
		{
			g_keys[wParam] = TRUE;					// If So, Mark It As TRUE
			return 0;								// Jump Back
		}

		case WM_KEYUP:								// Has A Key Been Released?
		{
			g_keys[wParam] = FALSE;					// If So, Mark It As FALSE
			return 0;								// Jump Back
		}

		case WM_SIZE:								// Resize The OpenGL Window
		{
			ReSizeGLScene(LOWORD(lParam),HIWORD(lParam));  // LoWord=Width, HiWord=Height
			return 0;								// Jump Back
		}

	}

	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

int WINAPI WinMain(	HINSTANCE	hInstance,			// Instance
					HINSTANCE	hPrevInstance,		// Previous Instance
					LPSTR		lpCmdLine,			// Command Line Parameters
					int			nCmdShow)			// Window Show State
{
	MSG		msg;									// Windows Message Structure

	// Create Our OpenGL Window
	if (!CreateGLWindow("Protein Channel Viewer",windowWidth,windowHeight,16,fullscreen))
	{
		return 0;									// Quit If Window Was Not Created
	}

	//Open a console window
    AllocConsole();
    //Connect console output
    HANDLE handle_out = GetStdHandle(STD_OUTPUT_HANDLE);
    int hCrt          = _open_osfhandle((long) handle_out, _O_TEXT);
    FILE* hf_out      = _fdopen(hCrt, "w");
    setvbuf(hf_out, NULL, _IONBF, 1);
    *stdout = *hf_out;
    //Connect console input
    HANDLE handle_in = GetStdHandle(STD_INPUT_HANDLE);
    hCrt             = _open_osfhandle((long) handle_in, _O_TEXT);
    FILE* hf_in      = _fdopen(hCrt, "r");
    setvbuf(hf_in, NULL, _IONBF, 128);
    *stdin = *hf_in;
    //Set console title
    SetConsoleTitle("PCV - Console");
    //Move console window
    hWnd_console = GetConsoleWindow();
    MoveWindow(hWnd_console,0,515,680,510,TRUE);

	while(!g_exit)									// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))	// Is There A Message Waiting?
		{
			if (msg.message==WM_QUIT)				// Have We Received A Quit Message?
			{
				g_exit=TRUE;							// If So done=TRUE
			}
			else									// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);				// Translate The Message
				DispatchMessage(&msg);				// Dispatch The Message
			}
		}
		else										// If There Are No Messages
		{
			// Draw The Scene.  Watch For ESC Key And Quit Messages From ProgramCycle()
			if (active)								// Program Active?
			{
				if (g_keys[VK_ESCAPE])				// Was ESC Pressed?
				{
					g_state=stEXIT;					// ESC Signalled A Quit
				}
				else                                // Not Time To Quit, Update Screen
				{
                    iLastTickCount=iTickCount;      // Save Last Tick Count
					iTickCount=GetTickCount();      // Get Current Tick Count
					if (iTickCount!=iLastTickCount) FPS=(1000/(iTickCount-iLastTickCount)); // Calculate FPS
					g_cycleTime=iTickCount-iLastTickCount; // Get the Time of Previous Cycle

                    if((g_keys[67]||g_keys[82])&&fullscreen)      // If console needs input and fullscreen is enabled, turn off fullscreen
                    {
                        fullscreen=false;
                        KillGLWindow();
                        windowWidth=640;
                        windowHeight=480;
                        // Recreate Our OpenGL Window
                        if (!CreateGLWindow("Protein Channel Viewer",windowWidth,windowHeight,16,fullscreen))
                        {
                            return 0;						// Quit If Window Was Not Created
                        }
                        pProtein->setWindowVal(windowWidth,windowHeight,screenWidth,screenHeight,hWnd_console,hWnd); // Tell Window Size to Protein
                    }

                    keyCheck();                     // Get User Input
					ProgramCycle();					// Updating
					SwapBuffers(hDC);				// Swap Buffers (Double Buffering)
				}
			}

			if (g_keys[VK_F1])						// Is F1 Being Pressed?
			{
				g_keys[VK_F1]=FALSE;				// If So Make Key FALSE
				KillGLWindow();						// Kill Our Current Window
				if(!fullscreen)
				{
                    fullscreen=TRUE;                // Toogle Fullscreen
                    windowWidth=GetSystemMetrics(SM_CXSCREEN); // Get Width of Screen
                    windowHeight=GetSystemMetrics(SM_CYSCREEN);// Get Height of Screen
				}
				else
				{
                    fullscreen=FALSE;               // Toogle Fullscreen
                    windowWidth=640;
                    windowHeight=480;
				}

				// Recreate Our OpenGL Window
				if (!CreateGLWindow("Protein Channel Viewer",windowWidth,windowHeight,16,fullscreen))
				{
					return 0;						// Quit If Window Was Not Created
				}
				pProtein->setWindowVal(windowWidth,windowHeight,screenWidth,screenHeight,hWnd_console,hWnd); // Tell Window Size to Protein
			}
		}
	}

	// Shutdown
	KillGLWindow();									// Kill The Window
	return (msg.wParam);							// Exit The Program
}
