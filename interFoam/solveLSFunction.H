//NOTE: SOLVELSFUNCTION
// 重初始化 LS方程

Info << "解重初始化方程" << nl << endl;

psi == psi0;

// corr : 循环次数
for(int corr = 0; corr < int(epsilon.value()/deltaTau.value()); corr++)
{
    psi = psi + psi0 / mag(psi0) * (double(1) - mag(fvc::grad(psi) *
    dimChange)) * deltaTau;
    psi.correctBoundaryConditions();
}

// 更新狄拉克函数
forAll(mesh.cells(), celli)
{
    if(mag(psi[celli]) > epsilon.value())
        delta[celli] = 0;
    else
        delta[celli] = 1.0 / (2.0 * epsilon.value()) *
        (1.0 + Foam::cos(M_PI * psi[celli] / epsilon.value()));
};

// 更新阶跃函数
forAll(mesh.cells(), celli)
{
    if(psi[celli] < -epsilon.value())
        H[celli] = double(0);
    else if(epsilon.value() < psi[celli])
        H[celli] = double(1);
    else
        H[celli] = double(1.0)/double(2.0)*(double(1.0)+psi[celli]/
        epsilon.value()+Foam::sin(M_PI*psi[celli]/epsilon.value())/M_PI);
};

// alpha1 = H
// alpha2 = 1-alpha1
// H限定在0~1
const volScalarField limitedH
(
    "limitedH",
    min(max(H, scalar(0)), scalar(1))
);

// 更新密度
rho == limitedH * rho1 + (1 - limitedH) * rho2;

const_cast<volScalarField&>(mixture.nu()()) = limitedH * nu1 + (1.0 - limitedH) * nu2;

H == limitedH;
