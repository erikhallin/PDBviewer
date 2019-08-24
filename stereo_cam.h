#ifndef STEREO_CAM_H
#define STEREO_CAM_H

#include <math.h>
#include <gl\gl.h>
#define _piover180 0.0174532925

class stereo_cam
{
    public:
        stereo_cam(float Convergence,
                     float EyeSeparation,
                     float AspectRatio,
                     float FOV,
                     float NearClippingDistance,
                     float FarClippingDistance);


        void ApplyLeftFrustum(void);
        void ApplyRightFrustum(void);

    private:
        float mConvergence;
        float mEyeSeparation;
        float mAspectRatio;
        float mFOV;
        float mNearClippingDistance;
        float mFarClippingDistance;
};

#endif // STEREO_CAM_H
