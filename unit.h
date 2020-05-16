// Unit tests for the math expression parser
// @author Stepan Belousov
#include <cmath>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "expression.h"

using namespace std;

class ExpressionTest: public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(ExpressionTest);
    CPPUNIT_TEST(testUnaryOperators);
    CPPUNIT_TEST(testBinaryOperators);
    CPPUNIT_TEST(testFunctions);
    CPPUNIT_TEST(testComplexExpr);
    CPPUNIT_TEST(testVariables);
    CPPUNIT_TEST(testWrongFormat);
    CPPUNIT_TEST(testBadResult);
    CPPUNIT_TEST(testUnknownVariable);
    CPPUNIT_TEST_SUITE_END();

    string varsAB, varsABC, varsAXT, varsX, varsConst, varsA, varsB, varsABCD;
    string varsBad1, varsBad2, varsBad3;

    static const double maxAbsError = 1e-5;
public:
    void setUp() {
        varsAB = "a = 1 b = 2";
        varsABC = "a = 1 b = 2 c = 3";
        varsAXT = "a = 1 t         = -0.1 x = 0.25";
        varsX = "x = 1.000";
        varsConst = "const = 1 const1 = 2 const2 = 3";
        varsA = "a = 1 aa = 2";
        varsB = "b = 2";
        varsABCD = "a = 1 b = 2 c = 3 d = 4";
        varsBad1 = "1 = 2";
        varsBad2 = "a ~ 1";
        varsBad3 = "a = e0";
    }

    bool resCheck(double res, double right) {
        return fabs(res - right) < maxAbsError;
    }

    bool wrongFormatCheck(const string& s, const string& vars = "") {
        try {
            Expression ex(s);
            ex.eval(vars);
        } catch (const ExpressionWrongFormat&) {
            return true;
        } catch (...) {}
        return false;
    }

    bool badResultCheck(const string& s, const string& vars = "") {
        try {
            Expression ex(s);
            ex.eval(vars);
        } catch (const ExpressionBadResult&) {
            return true;
        } catch (...) {}
        return false;
    }

    bool unknownVariableCheck(const string& s, const string& vars = "") {
        try {
            Expression ex(s);
            ex.eval(vars);
        } catch (const ExpressionUnknownVariable&) {
            return true;
        } catch (...) {}
        return false;
    }

    void testUnaryOperators() {
        CPPUNIT_ASSERT(resCheck(Expression("-(a+b)").eval(varsAB), -3.0));
        CPPUNIT_ASSERT(resCheck(Expression("-(-(-(-(-x))))").eval(varsX), -1.0));
        CPPUNIT_ASSERT(resCheck(Expression("-5.55").eval(), -5.55));
        CPPUNIT_ASSERT(resCheck(Expression("(+a)+(+b)").eval(varsAB), 3.0));
        CPPUNIT_ASSERT(resCheck(Expression("+(-(+(-(+x))))").eval(varsX), 1.0));
        CPPUNIT_ASSERT(resCheck(Expression("1-(+2)").eval(), -1.0));
        CPPUNIT_ASSERT(resCheck(Expression("+1").eval(), +1));
        CPPUNIT_ASSERT(resCheck(Expression("-(+1)").eval(), -1));
        CPPUNIT_ASSERT(resCheck(Expression("-(+1)*2").eval(),-2));
        CPPUNIT_ASSERT(resCheck(Expression("-(+2)*sqrt(4)").eval(), -4));
        CPPUNIT_ASSERT(resCheck(Expression("3-(+x)").eval(varsX), 2));
        CPPUNIT_ASSERT(resCheck(Expression("+1*3").eval(), 3));
        CPPUNIT_ASSERT(resCheck(Expression("-1").eval(), -1));
        CPPUNIT_ASSERT(resCheck(Expression("-(-1)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("-(-1)*2").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("-(-2)*sqrt(4)").eval(), 4));
        CPPUNIT_ASSERT(resCheck(Expression("-pi").eval(), -3.141592653589793));
        CPPUNIT_ASSERT(resCheck(Expression("-x").eval(varsX), -1));
        CPPUNIT_ASSERT(resCheck(Expression("-(x)").eval(varsX), -1));
        CPPUNIT_ASSERT(resCheck(Expression("-(-x)").eval(varsX), 1));
        CPPUNIT_ASSERT(resCheck(Expression("-(-x)*2").eval(varsX), 2));
        CPPUNIT_ASSERT(resCheck(Expression("-(8)").eval(), -8));
        CPPUNIT_ASSERT(resCheck(Expression("-8").eval(), -8));
        CPPUNIT_ASSERT(resCheck(Expression("-(2+1)").eval(), -3));
        CPPUNIT_ASSERT(resCheck(Expression("-sin(8)").eval(), -0.9893582466233818));
        CPPUNIT_ASSERT(resCheck(Expression("3-(-x)").eval(varsX), 4));
        CPPUNIT_ASSERT(resCheck(Expression("-1*3").eval(), -3));
        CPPUNIT_ASSERT(resCheck(Expression("(-3)^2").eval(), 9));
        CPPUNIT_ASSERT(resCheck(Expression("-(-2^2)").eval(), 4));
        CPPUNIT_ASSERT(resCheck(Expression("3+(-3^2)").eval(), -6));
    }

    void testBinaryOperators() {
        CPPUNIT_ASSERT(resCheck(Expression("a + b").eval(varsAB), 3.0));
        CPPUNIT_ASSERT(resCheck(Expression("x + 2").eval(varsX), 3.0));
        CPPUNIT_ASSERT(resCheck(Expression("5 + 6").eval(), 11.0));
        CPPUNIT_ASSERT(resCheck(Expression("a - b").eval(varsAB), -1.0));
        CPPUNIT_ASSERT(resCheck(Expression("x - 1").eval(varsX), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("0 - 3").eval(), -3.0));
        CPPUNIT_ASSERT(resCheck(Expression("a * b").eval(varsAB), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("x * 0.0").eval(varsX), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("-0.5 * 0.25").eval(), -0.125));
        CPPUNIT_ASSERT(resCheck(Expression("a / b").eval(varsAB), 0.5));
        CPPUNIT_ASSERT(resCheck(Expression("x / 1.0").eval(varsX), 1.0));
        CPPUNIT_ASSERT(resCheck(Expression("-1.0 / 0.5").eval(), -2.0));
        CPPUNIT_ASSERT(resCheck(Expression("1 / 2 / 3").eval(), 1.0 / 6.0));
        CPPUNIT_ASSERT(resCheck(Expression("b ^ a ^ b").eval(varsAB), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("4.5 ^ x").eval(varsX), 4.5));
        CPPUNIT_ASSERT(resCheck(Expression("100 ^ 0").eval(), 1.0));
        CPPUNIT_ASSERT(resCheck(Expression("2 ^ 2 ^ 3").eval(), 256.0));
    }

    void testFunctions() {
        CPPUNIT_ASSERT(resCheck(Expression("sin(x)").eval(varsX), 0.8414709848078965));
        CPPUNIT_ASSERT(resCheck(Expression("cos(x)").eval(varsX), 0.5403023058681398));
        CPPUNIT_ASSERT(resCheck(Expression("tg(x)").eval(varsX), 1.5574077246549023));
        CPPUNIT_ASSERT(resCheck(Expression("ctg(x)").eval(varsX), 0.6420926159343306));
        CPPUNIT_ASSERT(resCheck(Expression("arcsin(x)").eval(varsX), 1.5707963267948966));
        CPPUNIT_ASSERT(resCheck(Expression("arccos(-x)").eval(varsX), 3.141592653589793));
        CPPUNIT_ASSERT(resCheck(Expression("arctg(x)").eval(varsX), 0.7853981633974483));
        CPPUNIT_ASSERT(resCheck(Expression("arcctg(x*0)").eval(varsX), 1.5707963267948966));
        CPPUNIT_ASSERT(resCheck(Expression("exp(x)").eval(varsX), 2.718281828459045));
        CPPUNIT_ASSERT(resCheck(Expression("ln(x+1)").eval(varsX), 0.6931471805599453));
        CPPUNIT_ASSERT(resCheck(Expression("log10(x+11)").eval(varsX), 1.0791812460476247));
        CPPUNIT_ASSERT(resCheck(Expression("abs(-x)").eval(varsX), 1.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt(x+1)").eval(varsX), 1.4142135623730951));
        CPPUNIT_ASSERT(resCheck(Expression("sinh(x)").eval(varsX), 1.1752011936438014));
        CPPUNIT_ASSERT(resCheck(Expression("cosh(x)").eval(varsX), 1.5430806348152437));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt(4)").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt(4.)").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("pi").eval(), 3.14159265));
        CPPUNIT_ASSERT(resCheck(Expression("sin(pi)").eval(), 0));
        CPPUNIT_ASSERT(resCheck(Expression("sin(pi/2)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("sin(0)").eval(), 0));
        CPPUNIT_ASSERT(resCheck(Expression("sin(sin(0.1))").eval(), 0.09966766));
        CPPUNIT_ASSERT(resCheck(Expression("cos(0)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("exp(0)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("exp(3)").eval(), 20.085537));
        CPPUNIT_ASSERT(resCheck(Expression("ln(exp(2))").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("log10(100)").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("log2(8)").eval(), 3));
        CPPUNIT_ASSERT(resCheck(Expression("sinh(0)").eval(), 0));
        CPPUNIT_ASSERT(resCheck(Expression("sinh(1)").eval(), 1.1752012));
        CPPUNIT_ASSERT(resCheck(Expression("sinh(sinh(1.1))").eval(), 1.76973464));
        CPPUNIT_ASSERT(resCheck(Expression("cosh(0)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("cosh(1)").eval(), 1.5430806));
        CPPUNIT_ASSERT(resCheck(Expression("arcsin(0.5)").eval(), 0.523599));
        CPPUNIT_ASSERT(resCheck(Expression("arccos(0.5)").eval(), 1.0471976));
        CPPUNIT_ASSERT(resCheck(Expression("tg(0.5)").eval(), 0.54630249));
        CPPUNIT_ASSERT(resCheck(Expression("arctg(0.5)").eval(), 0.46364761));
    }

    void testComplexExpr() {
        CPPUNIT_ASSERT(resCheck(Expression("sin(2*pi*x)*exp(-4*pi^2*a*t)").eval(varsAXT), 51.82339874179074));
        CPPUNIT_ASSERT(resCheck(Expression("0.5 * (1. + 0.2 - 1e-1)").eval(), 0.55));
        CPPUNIT_ASSERT(resCheck(Expression("log10(3) - ln(3) / ln(10)").eval(), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("sin(cos(x)) + cos(sin(x))").eval(varsX), 1.1807620039164297));
        CPPUNIT_ASSERT(resCheck(Expression("a + b * c").eval(varsABC), 7.0));
        CPPUNIT_ASSERT(resCheck(Expression("a * b + c").eval(varsABC), 5.0));
        CPPUNIT_ASSERT(resCheck(Expression("x + 2 * 3").eval(varsX), 7));
        CPPUNIT_ASSERT(resCheck(Expression("3 + x + 2").eval(varsX), 6));
        CPPUNIT_ASSERT(resCheck(Expression("x * 2 * 3").eval(varsX), 6));
        CPPUNIT_ASSERT(resCheck(Expression("x + (2 + 3)").eval(varsX), 6));
        CPPUNIT_ASSERT(resCheck(Expression("x * (2 * 3)").eval(varsX), 6));
        CPPUNIT_ASSERT(resCheck(Expression("x + 2 ^ 3").eval(varsX), 9));
        CPPUNIT_ASSERT(resCheck(Expression("x ^ (-2)").eval(varsX), 1));
        CPPUNIT_ASSERT(resCheck(Expression("x * 2 * (3 * 4)").eval(varsX), 24));
        CPPUNIT_ASSERT(resCheck(Expression("x + 2 + (3 + 4)").eval(varsX), 10));
        CPPUNIT_ASSERT(resCheck(Expression("a ^ a - b ^ b").eval(varsAB), -3));
        CPPUNIT_ASSERT(resCheck(Expression("1+2*3").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("1+2*(-3)").eval(), -5));
        CPPUNIT_ASSERT(resCheck(Expression("1-2*3").eval(), -5));
        CPPUNIT_ASSERT(resCheck(Expression("1+4/2").eval(), 3));
        CPPUNIT_ASSERT(resCheck(Expression("1-4/2").eval(), -1));
        CPPUNIT_ASSERT(resCheck(Expression("3*2^4").eval(), 48));
        CPPUNIT_ASSERT(resCheck(Expression("4*2^(-2)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("3-2^4").eval(), -13));
        CPPUNIT_ASSERT(resCheck(Expression("3+2^4").eval(), 19));
        CPPUNIT_ASSERT(resCheck(Expression("16/2^3").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("2 * 2 ^4").eval(), 32));
        CPPUNIT_ASSERT(resCheck(Expression("(1)").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("(1+2)*3").eval(), 9));
        CPPUNIT_ASSERT(resCheck(Expression("(1-2)*3").eval(), -3));
        CPPUNIT_ASSERT(resCheck(Expression("(1+2)^2").eval(), 9));
        CPPUNIT_ASSERT(resCheck(Expression("(2)^(1+2)").eval(), 8));
        CPPUNIT_ASSERT(resCheck(Expression("((1))").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("0").eval(), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("0.0").eval(), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("1-1").eval(), 0.0));
        CPPUNIT_ASSERT(resCheck(Expression("1e1").eval(), 10));
        CPPUNIT_ASSERT(resCheck(Expression("1E1").eval(), 10));
        CPPUNIT_ASSERT(resCheck(Expression("1.e1").eval(), 10));
        CPPUNIT_ASSERT(resCheck(Expression("1.E1").eval(), 10));
        CPPUNIT_ASSERT(resCheck(Expression("3+(-3^2)").eval(), -6.0));
        CPPUNIT_ASSERT(resCheck(Expression("-2^2").eval(), -4.0));
        CPPUNIT_ASSERT(resCheck(Expression("(1+ 2*x)").eval(varsX), 3.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt((4))").eval(), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt((2)+2)").eval(), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt(2+(2))").eval(), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt(x+(3))").eval(varsX), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("sqrt((3)+x)").eval(varsX), 2.0));
        CPPUNIT_ASSERT(resCheck(Expression("2*b*5").eval(varsB), 20));
        CPPUNIT_ASSERT(resCheck(Expression("2*b*5 + 4*b").eval(varsB), 28));
        CPPUNIT_ASSERT(resCheck(Expression("2*x/3").eval(varsX), 2.0/3.0));
        CPPUNIT_ASSERT(resCheck(Expression("3+b").eval(varsB), 5));
        CPPUNIT_ASSERT(resCheck(Expression("b+3").eval(varsB), 5));
        CPPUNIT_ASSERT(resCheck(Expression("b*3+2").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("3*b+2").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("2+b*3").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("2+3*b").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("b+3*b").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("3*b+b").eval(varsB), 8));
        CPPUNIT_ASSERT(resCheck(Expression("2+b*3+b").eval(varsB), 10));
        CPPUNIT_ASSERT(resCheck(Expression("b+2+b*3").eval(varsB), 10));
        CPPUNIT_ASSERT(resCheck(Expression("(2*b+1)*4").eval(varsB), 20));
        CPPUNIT_ASSERT(resCheck(Expression("4*(2*b+1)").eval(varsB), 20));
        CPPUNIT_ASSERT(resCheck(Expression("1+2-3*4/5^6").eval(), 2.99923));
        CPPUNIT_ASSERT(resCheck(Expression("1^2/3*4-5+6").eval(), 2.33333333));
        CPPUNIT_ASSERT(resCheck(Expression("1+2*3").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("1+2*3").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("(1+2)*3").eval(), 9));
        CPPUNIT_ASSERT(resCheck(Expression("(1+2)*(-3)").eval(), -9));
        CPPUNIT_ASSERT(resCheck(Expression("2/4").eval(), 0.5));
        CPPUNIT_ASSERT(resCheck(Expression("1+2").eval(), 3));
        CPPUNIT_ASSERT(resCheck(Expression("1+2+3").eval(), 6));
        CPPUNIT_ASSERT(resCheck(Expression("2*3").eval(), 6));
        CPPUNIT_ASSERT(resCheck(Expression("2*3*4").eval(), 24));
        CPPUNIT_ASSERT(resCheck(Expression("3-2").eval(), 1));
        CPPUNIT_ASSERT(resCheck(Expression("3-2-1").eval(), 0));
        CPPUNIT_ASSERT(resCheck(Expression("4/2").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("16/4/2").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("4*3/2").eval(), 6));
        CPPUNIT_ASSERT(resCheck(Expression("4/2*3").eval(), 6));
        CPPUNIT_ASSERT(resCheck(Expression("1+2-3").eval(), 0));
        CPPUNIT_ASSERT(resCheck(Expression("1-2+3").eval(), 2));
        CPPUNIT_ASSERT(resCheck(Expression("exp(ln(7))").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("exp(1)^ln(7)").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("exp(1)^(ln(7))").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("(exp(1)^(ln(7)))").eval(), 7));
        CPPUNIT_ASSERT(resCheck(Expression("1-(exp(1)^(ln(7)))").eval(), -6));
        CPPUNIT_ASSERT(resCheck(Expression("2*(exp(1)^(ln(7)))").eval(), 14));
        CPPUNIT_ASSERT(resCheck(Expression("10^ln(5)").eval(), pow(10.0, log(5.0))));
        CPPUNIT_ASSERT(resCheck(Expression("10^log10(5)").eval(), 5));
        CPPUNIT_ASSERT(resCheck(Expression("2^log2(4)").eval(), 4));
        CPPUNIT_ASSERT(resCheck(Expression("-(sin(0)+1)").eval(), -1));
        CPPUNIT_ASSERT(resCheck(Expression("-(2^1.1)").eval(), -2.14354692));
        CPPUNIT_ASSERT(resCheck(Expression("(cos(2.41)/b)").eval(varsB), -0.372056));
        CPPUNIT_ASSERT(resCheck(Expression("(1*(2*(3*(4*(5*(6*(a+b)))))))").eval(varsAB), 2160));
        CPPUNIT_ASSERT(resCheck(Expression("(1*(2*(3*(4*(5*(6*(7*(a+b))))))))").eval(varsAB), 15120));
        CPPUNIT_ASSERT(resCheck(Expression("(a/((((b+(((exp(1)*(((((pi*((((3.45*((pi+a)+pi))+b)+b)*a))+0.68)+exp(1))+a)/a))+a)+b))+b)*a)-pi))").eval(varsAB), 0.00377999));
        CPPUNIT_ASSERT(resCheck(Expression("1+2-3*4/5^6*(2*(1-5+(3*7^9)*(4+6*7-3)))+12").eval(), -7995810.09926));
        CPPUNIT_ASSERT(resCheck(Expression("(arctg(sin((((((((((((((((pi/cos((a/((((0.53-b)-pi)*exp(1))/b))))+2.51)+a)-0.54)/0.98)+b)*b)+exp(1))/a)+b)+a)+b)+pi)/exp(1))+a)))*2.77)").eval(varsAB), -2.16995656));
        CPPUNIT_ASSERT(resCheck(Expression(
        "(((-9))-exp(1)/(((((((pi-(((-7)+(-3)/4/exp(1)))))/(((-5))-2)-((pi+(-0))*(sqrt((exp(1)+exp(1)))*(-8))*(((-pi)+(-pi)-(-9)*(6*5))"
        "/(-exp(1))-exp(1)))/2)/((((sqrt(2/(-exp(1))+6)-(4-2))+((5/(-2))/(1*(-pi)+3))/8)*pi*((pi/((-2)/(-6)*1*(-1))*(-6)+(-exp(1))))))/"
        "((exp(1)+(-2)+(-exp(1))*((((-3)*9+(-exp(1))))+(-9)))))))-((((exp(1)-7+(((5/pi-(3/1+pi)))))/exp(1))/(-5))/(sqrt((((((1+(-7))))+((((-"
        "exp(1))*(-exp(1))))-8))*(-5)/((-exp(1))))*(-6)-((((((-2)-(-9)-(-exp(1))-1)/3))))/(sqrt((8+(exp(1)-((-6))+(9*(-9))))*(((3+2-8))*(7+6"
        "+(-5))+((0/(-exp(1))*(-pi))+7)))+(((((-exp(1))/exp(1)/exp(1))+((-6)*5)*exp(1)+(3+(-5)/pi))))+pi))/sqrt((((9))+((((pi))-8+2))+pi))/exp(1)"
        "*4)*((-5)/(((-pi))*(sqrt(exp(1))))))-(((((((-exp(1))*(exp(1))-pi))/4+(pi)*(-9)))))))+(-pi)").eval(), -12.23016549));
    }

    void testVariables() {
        CPPUNIT_ASSERT(resCheck(Expression("const").eval(varsConst), 1));
        CPPUNIT_ASSERT(resCheck(Expression("const1").eval(varsConst), 2));
        CPPUNIT_ASSERT(resCheck(Expression("const2").eval(varsConst), 3));
        CPPUNIT_ASSERT(resCheck(Expression("2*const").eval(varsConst), 2));
        CPPUNIT_ASSERT(resCheck(Expression("2*const1").eval(varsConst), 4));
        CPPUNIT_ASSERT(resCheck(Expression("2*const2").eval(varsConst), 6));
        CPPUNIT_ASSERT(resCheck(Expression("2*const+1").eval(varsConst), 3));
        CPPUNIT_ASSERT(resCheck(Expression("2*const1+1").eval(varsConst), 5));
        CPPUNIT_ASSERT(resCheck(Expression("2*const2+1").eval(varsConst), 7));
        CPPUNIT_ASSERT(resCheck(Expression("a").eval(varsA), 1));
        CPPUNIT_ASSERT(resCheck(Expression("aa").eval(varsA), 2));
        CPPUNIT_ASSERT(resCheck(Expression("2*a").eval(varsA), 2));
        CPPUNIT_ASSERT(resCheck(Expression("2*aa").eval(varsA), 4));
        CPPUNIT_ASSERT(resCheck(Expression("2*a-1").eval(varsA), 1));
        CPPUNIT_ASSERT(resCheck(Expression("2*aa-1").eval(varsA), 3));
        CPPUNIT_ASSERT(resCheck(Expression("a + b + c + d").eval(varsABCD), 10));
    }

    void testWrongFormat() {
        CPPUNIT_ASSERT(wrongFormatCheck(""));
        CPPUNIT_ASSERT(wrongFormatCheck("+"));
        CPPUNIT_ASSERT(wrongFormatCheck("-"));
        CPPUNIT_ASSERT(wrongFormatCheck("*"));
        CPPUNIT_ASSERT(wrongFormatCheck("/"));
        CPPUNIT_ASSERT(wrongFormatCheck("^"));
        CPPUNIT_ASSERT(wrongFormatCheck("("));
        CPPUNIT_ASSERT(wrongFormatCheck(")"));
        CPPUNIT_ASSERT(wrongFormatCheck("()"));
        CPPUNIT_ASSERT(wrongFormatCheck("((1)"));
        CPPUNIT_ASSERT(wrongFormatCheck("(1))"));
        CPPUNIT_ASSERT(wrongFormatCheck("(1"));
        CPPUNIT_ASSERT(wrongFormatCheck("1)"));
        CPPUNIT_ASSERT(wrongFormatCheck("1)))((("));
        CPPUNIT_ASSERT(wrongFormatCheck(")()(1"));
        CPPUNIT_ASSERT(wrongFormatCheck("0a"));
        CPPUNIT_ASSERT(wrongFormatCheck("9a"));
        CPPUNIT_ASSERT(wrongFormatCheck("123abc"));
        CPPUNIT_ASSERT(wrongFormatCheck("0e"));
        CPPUNIT_ASSERT(wrongFormatCheck("0x"));
        CPPUNIT_ASSERT(wrongFormatCheck("0x11"));
        CPPUNIT_ASSERT(wrongFormatCheck("0exp"));
        CPPUNIT_ASSERT(wrongFormatCheck("0 exp"));
        CPPUNIT_ASSERT(wrongFormatCheck("0cos"));
        CPPUNIT_ASSERT(wrongFormatCheck("0 cos"));        
        CPPUNIT_ASSERT(wrongFormatCheck("sin 0"));
        CPPUNIT_ASSERT(wrongFormatCheck("1 2 +"));
        CPPUNIT_ASSERT(wrongFormatCheck("+ 1 2"));
        CPPUNIT_ASSERT(wrongFormatCheck("0 -"));
        CPPUNIT_ASSERT(wrongFormatCheck("0 +"));
        CPPUNIT_ASSERT(wrongFormatCheck("++0"));
        CPPUNIT_ASSERT(wrongFormatCheck("--0"));
        CPPUNIT_ASSERT(wrongFormatCheck("++x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("--x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("/1"));
        CPPUNIT_ASSERT(wrongFormatCheck("*1"));
        CPPUNIT_ASSERT(wrongFormatCheck("2^-2"));
        CPPUNIT_ASSERT(wrongFormatCheck("2^+2"));
        CPPUNIT_ASSERT(wrongFormatCheck("1--1"));
        CPPUNIT_ASSERT(wrongFormatCheck("1++1"));
        CPPUNIT_ASSERT(wrongFormatCheck("1**1"));
        CPPUNIT_ASSERT(wrongFormatCheck("1//1"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin^2(x)", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("sin-(x)", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("sin x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("1+"));
        CPPUNIT_ASSERT(wrongFormatCheck("1-"));
        CPPUNIT_ASSERT(wrongFormatCheck("x+", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("x-", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("x*", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("x/", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("?x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("!x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("x?", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("x!", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("1,"));
        CPPUNIT_ASSERT(wrongFormatCheck("a,"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(8),"));
        CPPUNIT_ASSERT(wrongFormatCheck("(sin(8)),"));
        CPPUNIT_ASSERT(wrongFormatCheck("{2}"));
        CPPUNIT_ASSERT(wrongFormatCheck("{ 2 }"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin{2}"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin{2}"));
        CPPUNIT_ASSERT(wrongFormatCheck("(sin){2}"));
        CPPUNIT_ASSERT(wrongFormatCheck("(2+"));
        CPPUNIT_ASSERT(wrongFormatCheck("2++4"));
        CPPUNIT_ASSERT(wrongFormatCheck("2+-4"));
        CPPUNIT_ASSERT(wrongFormatCheck("(2+)"));
        CPPUNIT_ASSERT(wrongFormatCheck("--2"));
        CPPUNIT_ASSERT(wrongFormatCheck("()"));
        CPPUNIT_ASSERT(wrongFormatCheck("5+()"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(cos)"));
        CPPUNIT_ASSERT(wrongFormatCheck("5t6"));
        CPPUNIT_ASSERT(wrongFormatCheck("5 t 6"));
        CPPUNIT_ASSERT(wrongFormatCheck("8*"));
        CPPUNIT_ASSERT(wrongFormatCheck(",3"));
        CPPUNIT_ASSERT(wrongFormatCheck("3,5"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(8,8)"));
        CPPUNIT_ASSERT(wrongFormatCheck("(7,8)"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin)"));
        CPPUNIT_ASSERT(wrongFormatCheck("a)"));
        CPPUNIT_ASSERT(wrongFormatCheck("pi)"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(())"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin()"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(1,1)"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(...)"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin."));
        CPPUNIT_ASSERT(wrongFormatCheck("(sin)(1)"));
        CPPUNIT_ASSERT(wrongFormatCheck("3+"));
        CPPUNIT_ASSERT(wrongFormatCheck("3+)"));
        CPPUNIT_ASSERT(wrongFormatCheck("3+()"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(3 4)"));
        CPPUNIT_ASSERT(wrongFormatCheck("(1+2"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(3)3"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(3)xyz"));
        CPPUNIT_ASSERT(wrongFormatCheck("sin(3)cos(3)"));
        CPPUNIT_ASSERT(wrongFormatCheck("a = 1"));
        CPPUNIT_ASSERT(wrongFormatCheck("a+b+c=10"));
        CPPUNIT_ASSERT(wrongFormatCheck("a=b=3"));        
        CPPUNIT_ASSERT(wrongFormatCheck("1 2"));
        CPPUNIT_ASSERT(wrongFormatCheck("(1) (2)"));
        CPPUNIT_ASSERT(wrongFormatCheck("1 + * 2"));
        CPPUNIT_ASSERT(wrongFormatCheck("&"));
        CPPUNIT_ASSERT(wrongFormatCheck("&x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("&1"));
        CPPUNIT_ASSERT(wrongFormatCheck("_"));
        CPPUNIT_ASSERT(wrongFormatCheck("_x", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("_1"));
        CPPUNIT_ASSERT(wrongFormatCheck("a * b 3"));
        CPPUNIT_ASSERT(wrongFormatCheck("1;2"));
        CPPUNIT_ASSERT(wrongFormatCheck("1; 2"));
        CPPUNIT_ASSERT(wrongFormatCheck("1 ; 2"));
        CPPUNIT_ASSERT(wrongFormatCheck(".1"));
        CPPUNIT_ASSERT(wrongFormatCheck("e1.0"));
        CPPUNIT_ASSERT(wrongFormatCheck("x 1.0", varsX));
        CPPUNIT_ASSERT(wrongFormatCheck("1", varsBad1));
        CPPUNIT_ASSERT(wrongFormatCheck("1", varsBad2));
        CPPUNIT_ASSERT(wrongFormatCheck("1", varsBad3));
    }

    void testBadResult() {
        CPPUNIT_ASSERT(badResultCheck("1/0"));
        CPPUNIT_ASSERT(badResultCheck("sqrt(-1)"));
        CPPUNIT_ASSERT(badResultCheck("ln(0)"));
        CPPUNIT_ASSERT(badResultCheck("log2(0)"));
        CPPUNIT_ASSERT(badResultCheck("log10(0)"));
        CPPUNIT_ASSERT(badResultCheck("ln(-1)"));
        CPPUNIT_ASSERT(badResultCheck("log2(-1)"));
        CPPUNIT_ASSERT(badResultCheck("log10(-1)"));
        CPPUNIT_ASSERT(badResultCheck("arcsin(2)"));
        CPPUNIT_ASSERT(badResultCheck("arccos(2)"));
    }

    void testUnknownVariable() {
        CPPUNIT_ASSERT(unknownVariableCheck("a"));
        CPPUNIT_ASSERT(unknownVariableCheck("aa"));
        CPPUNIT_ASSERT(unknownVariableCheck("a1"));
        CPPUNIT_ASSERT(unknownVariableCheck("a", varsX));
        CPPUNIT_ASSERT(unknownVariableCheck("a1", varsB));
        CPPUNIT_ASSERT(unknownVariableCheck("tan"));
        CPPUNIT_ASSERT(unknownVariableCheck("a+b", varsA));
        CPPUNIT_ASSERT(unknownVariableCheck("a+b+c+d+e", varsABCD));
        CPPUNIT_ASSERT(unknownVariableCheck("e"));
        CPPUNIT_ASSERT(unknownVariableCheck("1+x"));
        CPPUNIT_ASSERT(unknownVariableCheck("0*a"));
    }
};
