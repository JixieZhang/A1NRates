{
    gSystem->Load("libA1NRates.so");

    gSystem->AddIncludePath("-I$PWD/include");
    gInterpreter->AddIncludePath("$PWD/include");
    gSystem->AddIncludePath("-I$PWD/XSModel");
    gInterpreter->AddIncludePath("$PWD/XSModel");
}
