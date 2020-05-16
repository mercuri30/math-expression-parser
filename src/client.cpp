#include <iostream>
#include <string>
#include "expression.h"

using namespace std;

int main() {
    cout.precision(6);
    cout.setf(ios::fixed);

    string line;
    while(getline(cin, line)) {
        string expression, variables;
        size_t varSep = line.find('|');
        if (varSep != string::npos) {
            variables = line.substr(0, varSep);
            expression = line.substr(varSep + 1);
        } else {
            expression = line;
            variables = "";
        }

        try {
            Expression expr(expression);
            cout << expr.eval(variables) << endl;
        } catch (const std::exception& e) {
            cout << e.what() << endl;
        }
    }

    return 0;
}
