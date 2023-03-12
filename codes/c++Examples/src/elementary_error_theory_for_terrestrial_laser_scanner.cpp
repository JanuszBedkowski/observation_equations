#include <imgui.h>
#include <imgui_impl_glut.h>
#include <imgui_impl_opengl2.h>
#include <imgui_internal.h>

#include <GL/freeglut.h>
#include <Eigen/Eigen>

#include <elementary_error_theory_for_terrestrial_laser_scanner_jacobian.h>

#include <iostream>

static bool show_demo_window = true;
static bool show_another_window = false;
static ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
float rotate_x = 0.0, rotate_y = 0.0;
float translate_x, translate_y = 0.0;
float translate_z = -50.0;
const unsigned int window_width = 800;
const unsigned int window_height = 600;
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0; 
bool gui_mouse_down{ false };
float mouse_sensitivity = 1.0;

double sigma_r = 0.01;
double sigma_alpha = 0.01;
double sigma_theta = 0.01;
double measurement_r = 10.0;
double measurement_alpha = M_PI * 0.25;
double measurement_theta = M_PI * 0.25;

void my_display_code()
{
    // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
    if (show_demo_window)
        ImGui::ShowDemoWindow(&show_demo_window);

    // 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
    {
        static float f = 0.0f;
        static int counter = 0;

        ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

        ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
        ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
        ImGui::Checkbox("Another Window", &show_another_window);

        ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
        ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

        if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
            counter++;
        ImGui::SameLine();
        ImGui::Text("counter = %d", counter);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        ImGui::End();
    }

    // 3. Show another simple window.
    if (show_another_window)
    {
        ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
        ImGui::Text("Hello from another window!");
        if (ImGui::Button("Close Me"))
            show_another_window = false;
        ImGui::End();
    }
}

void reshape(int w, int h) {
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (GLfloat)w / (GLfloat)h, 0.01, 10000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void draw_ellipse(const Eigen::Matrix3d& covar, Eigen::Vector3d& mean, Eigen::Vector3f color, float nstd  = 3)
{
    Eigen::LLT<Eigen::Matrix<double,3,3> > cholSolver(covar);
    Eigen::Matrix3d transform = cholSolver.matrixL();

    const double pi = 3.141592;
    const double di = 0.02;
    const double dj = 0.04;
    const double du = di*2*pi;
    const double dv = dj*pi;
    glColor3f(color.x(), color.y(),color.z());

    for (double i = 0; i < 1.0; i+=di)  //horizonal
    {
        for (double j = 0; j < 1.0; j+=dj)  //vertical
        {
            double u = i*2*pi;      //0     to  2pi
            double v = (j-0.5)*pi;  //-pi/2 to pi/2

            const Eigen::Vector3d pp0( cos(v)* cos(u),cos(v) * sin(u),sin(v));
            const Eigen::Vector3d pp1(cos(v) * cos(u + du) ,cos(v) * sin(u + du) ,sin(v));
            const Eigen::Vector3d pp2(cos(v + dv)* cos(u + du) ,cos(v + dv)* sin(u + du) ,sin(v + dv));
            const Eigen::Vector3d pp3( cos(v + dv)* cos(u),cos(v + dv)* sin(u),sin(v + dv));
            Eigen::Vector3d tp0 = transform * (nstd*pp0) + mean;
            Eigen::Vector3d tp1 = transform * (nstd*pp1) + mean;
            Eigen::Vector3d tp2 = transform * (nstd*pp2) + mean;
            Eigen::Vector3d tp3 = transform * (nstd*pp3) + mean;

            glBegin(GL_LINE_LOOP);
            glVertex3dv(tp0.data());
            glVertex3dv(tp1.data());
            glVertex3dv(tp2.data());
            glVertex3dv(tp3.data());
            glEnd();
        }
    }
}

void display() {
    ImGuiIO& io = ImGui::GetIO();
    glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float ratio = float(io.DisplaySize.x) / float(io.DisplaySize.y);
    
    glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
        
    reshape((GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
    glTranslatef(translate_x, translate_y, translate_z);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 0.0, 1.0);
   
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);

    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
 
    
    double x = measurement_r * sin(measurement_alpha) * cos(measurement_theta);
    double y = measurement_r * sin(measurement_alpha) * sin(measurement_theta);
    double z = measurement_r * cos(measurement_theta);

    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(x, y, z);
    glEnd();
    
    Eigen::Matrix<double, 3, 3> j;

    elementary_error_theory_for_terrestrial_laser_scanner_jacobian(j, measurement_r, measurement_alpha, measurement_theta);
    
    Eigen::Matrix<double, 3, 3> cox_r_alpha_theta;
    cox_r_alpha_theta(0, 0) = sigma_r * sigma_r;
    cox_r_alpha_theta(0, 1) = 0.0;
    cox_r_alpha_theta(0, 2) = 0.0;

    cox_r_alpha_theta(1, 0) = 0.0;
    cox_r_alpha_theta(1, 1) = sigma_alpha * sigma_alpha;
    cox_r_alpha_theta(1, 2) = 0.0;

    cox_r_alpha_theta(1, 0) = 0.0;
    cox_r_alpha_theta(1, 1) = 0.0;
    cox_r_alpha_theta(1, 2) = sigma_theta * sigma_theta;

    
    Eigen::Matrix<double, 3, 3> cov_xyz = j * cox_r_alpha_theta * j.transpose();

    //std::cout << "---cov_xyz---" << std::endl;
    //std::cout << cov_xyz << std::endl;

    draw_ellipse(cov_xyz, Eigen::Vector3d(x, y, z), Eigen::Vector3f(0, 1, 0), 3);


    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGLUT_NewFrame();

    //my_display_code();
    if(ImGui::Begin("error propagation DEMO")){
        ImGui::InputDouble("sigma_r", &sigma_r, 0.0001, 0.001);
        ImGui::InputDouble("sigma_alpha", &sigma_alpha, 0.0001, 0.001);
        ImGui::InputDouble("sigma_theta", &sigma_theta, 0.0001, 0.001);
        ImGui::InputDouble("measurement_r", &measurement_r, 1.0, 1.0);
        ImGui::InputDouble("measurement_alpha", &measurement_alpha, 0.01, 0.1);
        ImGui::InputDouble("measurement_theta", &measurement_theta, 0.01, 0.1);
        ImGui::End();
    }


    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    
    glutSwapBuffers();
    glutPostRedisplay();
}

void motion(int x, int y) {
    ImGuiIO& io = ImGui::GetIO();
    io.MousePos = ImVec2((float)x, (float)y);

    if (!io.WantCaptureMouse)
    {
        float dx, dy;
        dx = (float)(x - mouse_old_x);
        dy = (float)(y - mouse_old_y);

        gui_mouse_down = mouse_buttons > 0;
        if (mouse_buttons & 1) {
            rotate_x += dy * 0.2f * mouse_sensitivity;
            rotate_y += dx * 0.2f * mouse_sensitivity;
        }
        if (mouse_buttons & 4) {
            translate_x += dx * 0.05f * mouse_sensitivity;
            translate_y -= dy * 0.05f * mouse_sensitivity;
        }
        mouse_old_x = x;
        mouse_old_y = y;
    }
    glutPostRedisplay();
}

void mouse(int glut_button, int state, int x, int y) {
    ImGuiIO& io = ImGui::GetIO();
    io.MousePos = ImVec2((float)x, (float)y);
    int button = -1;
    if (glut_button == GLUT_LEFT_BUTTON) button = 0;
    if (glut_button == GLUT_RIGHT_BUTTON) button = 1;
    if (glut_button == GLUT_MIDDLE_BUTTON) button = 2;
    if (button != -1 && state == GLUT_DOWN)
        io.MouseDown[button] = true;
    if (button != -1 && state == GLUT_UP)
        io.MouseDown[button] = false;

    if (!io.WantCaptureMouse)
    {
        if (state == GLUT_DOWN) {
            mouse_buttons |= 1 << glut_button;
        }
        else if (state == GLUT_UP) {
            mouse_buttons = 0;
        }
        mouse_old_x = x;
        mouse_old_y = y;
    }
}

void wheel(int button, int dir, int x, int y)
{
    if (dir > 0)
    {
        translate_z -= 0.05f * translate_z;
    }
    else
    {
        translate_z += 0.05f * translate_z;
    }
    return;
}

bool initGL(int* argc, char** argv) {
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("elementary_error_theory_for_terrestrial_laser_scanner");
    glutDisplayFunc(display);
    glutMotionFunc(motion);

    // default initialization
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    // viewport
    glViewport(0, 0, window_width, window_height);

    // projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (GLfloat)window_width / (GLfloat)window_height, 0.01, 10000.0);
    glutReshapeFunc(reshape);
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    ImGui::StyleColorsDark();
    ImGui_ImplGLUT_Init();
    ImGui_ImplGLUT_InstallFuncs();
    ImGui_ImplOpenGL2_Init();
    return true;
}

int main(int argc, char *argv[]){
    initGL(&argc, argv);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutMouseWheelFunc(wheel);
    glutMainLoop();

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
    return 0;
}