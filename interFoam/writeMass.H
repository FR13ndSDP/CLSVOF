{
    const scalarField &V = mesh.V();
    totalMass = gSum(rho * V) / totalMass0;

    // 如果文件不存在，首先创建
    if (massFilePtr.empty())
    {
        Info << "Creating mass file." << endl;

        // File update
        if (Pstream::master())
        {
            fileName totalMass;
            word name_ = "totalMass";

            // 打开新文件
            massFilePtr.reset(new OFstream("totalMass.dat"));
        }
    }

    if (Pstream::master())
    {
        massFilePtr() << runTime.timeName() << tab << totalMass << endl;
    }
}
