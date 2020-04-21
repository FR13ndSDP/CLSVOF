// 计算曲面法向量
volVectorField gradPsi(fvc::grad(psi));
surfaceVectorField gradPsif(fvc::interpolate(gradPsi));
surfaceVectorField nVecfv(gradPsif/(mag(gradPsif)+scalar(1.0e-6)/dimChange));
surfaceScalarField nVecf(nVecfv & mesh.Sf());

// 基于psi计算曲率
C == -fvc::div(nVecf);
