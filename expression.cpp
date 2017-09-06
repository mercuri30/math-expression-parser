#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <vector>
#include "expression.h"

using namespace std;

const string Expression::TokenTypeList[] = { "UNDEFINED", "PARENTHESIS", "UNDEFINED_OPERATOR", 
	"BINARY_OPERATOR", "UNARY_OPERATOR", "VARIABLE", "NUMBER", "FUNCTION" };

// List of supported functions, operators and constants
// If you wish to add smth here, don't forget to modify 
// applyFunction(), applyUnary() or applyBinary() as well
void Expression::initialize() {
	functions.insert("sin");
	functions.insert("cos");
	functions.insert("tg");
	functions.insert("ctg");
	functions.insert("arcsin");
	functions.insert("arccos");
	functions.insert("arctg");
	functions.insert("arcctg");
	functions.insert("exp");
	functions.insert("ln");
	functions.insert("log10");
	functions.insert("log2");
	functions.insert("abs");
	functions.insert("sqrt");
	functions.insert("sinh");
	functions.insert("cosh");
	operators.insert("+");
	operators.insert("-");
	operators.insert("*");
	operators.insert("/");
	operators.insert("^");
	variables["pi"] = pi;
}

double Expression::applyFunction(const string& s, double arg) const {
	if (s == "sin") return sin(arg);
	else if (s == "cos") return cos(arg);
	else if (s == "tg") return tan(arg);
	else if (s == "ctg") return 1.0 / tan(arg);
	else if (s == "arcsin") return asin(arg);
	else if (s == "arccos") return acos(arg);
	else if (s == "arctg") return atan(arg);
	else if (s == "arcctg") return pi / 2.0 - atan(arg);
	else if (s == "exp") return exp(arg);
	else if (s == "ln") return log(arg);
	else if (s == "log10") return log10(arg);
	else if (s == "log2") return log(arg) / log(2.0);
	else if (s == "abs") return fabs(arg);
	else if (s == "sqrt") return sqrt(arg);
	else if (s == "sinh") return sinh(arg);
	else if (s == "cosh") return cosh(arg);
	else throw ExpressionWrongFormat();
	return -1;
}

double Expression::applyUnary(const string& s, double arg) const {
	if (s == "+") return arg;
	else if (s == "-") return -arg;
	else throw ExpressionWrongFormat();
	return -1;
}

double Expression::applyBinary(const string& s, double arg1, double arg2) const {
	if (s == "+") return arg1 + arg2;
	else if (s == "-") return arg1 - arg2;
	else if (s == "*") return arg1 * arg2;
	else if (s == "/") return arg1 / arg2;
	else if (s == "^") return pow(arg1, arg2);
	else throw ExpressionWrongFormat();
	return -1;
}

// Parse variables and it's values from string and fill a map
// String must has a format like "x = 1 y = 2", etc.
void Expression::setVariables(const string& s) {
	istringstream iss(s);
	string s1, s2;
	double val;
	while (iss >> s1) {
		iss >> s2 >> val;
		if (!isVariable(s1) || iss.fail() || s2 != "=") {
			throw ExpressionWrongFormat();
		}
		variables[s1] = val;
	}
}

bool Expression::isOperator(const string& s) const {
	return (operators.find(s) != operators.end());
}

bool Expression::isDelimiter(char c) const {
	if (isOperator(string(1, c))) return true;
	return (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '(' || c == ')');
}

bool Expression::isVariable(const string& s) const {
	if (!isalpha(s[0])) return false;
	for (size_t i = 1;i < s.length();i++) {
		if (!isalpha(s[i]) && !isdigit(s[i])) return false;
	}
	return true;
}

Expression::TokenType Expression::getType(const Token& t) const {
	if (functions.find(t.s) != functions.end()) return FUNCTION;
	if (isVariable(t.s)) return VARIABLE;
	throw ExpressionWrongFormat();
	return UNDEFINED;
}

// Mark types for unary operators
void Expression::markUnaries(vector < Token > &v) const {
	if (v[0].type == UNDEFINED_OPERATOR) {
		v[0].type = UNARY_OPERATOR;
	}
	for (int i = static_cast<int>(v.size()) - 1;i >= 1;i--) {
		if (v[i].type != UNDEFINED_OPERATOR) continue;
		if (v[i - 1].type == UNDEFINED_OPERATOR || v[i - 1].type == FUNCTION) {
			throw ExpressionWrongFormat();
		}
		if (v[i - 1].s == "(") {
			v[i].type = UNARY_OPERATOR;
		} else {
			v[i].type = BINARY_OPERATOR;
		}
	}
}

// Split expression string into tokens
void Expression::split(const string& s, vector < Token > &result) {
	initialize();
	istringstream ss(s);
	string cur;
	while (true) {
		char c = static_cast<char>(ss.get());
		if (!ss) break;
		if (isDelimiter(c)) {
			if (cur != "") {
				result.push_back(Token(cur, getType(cur)));
				cur = "";
			}
			if (isOperator(string(1, c))) result.push_back(Token(string(1, c), UNDEFINED_OPERATOR));
			else if (c == '(' || c == ')') result.push_back(Token(string(1, c), PARENTHESIS));
		} else {
			if (cur == "" && isdigit(c)) {
				double num;
				ostringstream buffer;
				ss.unget();
				ss >> num;
				if (ss.fail()) {
					throw ExpressionWrongFormat();
				}
				buffer << num;
				result.push_back(Token(buffer.str(), NUMBER));
			} else {
				cur += c;
			}
		}
	}
	if (cur != "") {
		result.push_back(Token(cur, getType(cur)));
	}
}

// Get token priority for the Reverse Polish Notation
int Expression::priority(const Token& t) const {
	if (t.s == "(") return 0;
	if (t.type == FUNCTION) return 5;
	if (t.type == UNARY_OPERATOR) return 3;
	char c = t.s[0];
	switch (c) {
		case '+':
			return 1;
		case '-':
			return 1;
		case '*':
			return 2;
		case '/':
			return 2;
		case '^':
			return 4;
		default:
			throw ExpressionWrongFormat();
	}
	return -1;
}

// Convert expression to the Reverse Polish Notation
void Expression::convertExpression(vector < Token > &v, vector < Token > &result) const {
	markUnaries(v);
	stack < Token > st;
	int p;
	for (size_t i = 0;i < v.size();i++) {
		switch (v[i].type) {
			case PARENTHESIS:
				if (v[i].s == "(") {
					st.push(v[i]);
				} else {
					while (!st.empty() && st.top().s != "(") {
						result.push_back(st.top());
						st.pop();
					}
					if (st.empty()) {
						throw ExpressionWrongFormat();
					}
					st.pop();
				}
				break;
			case VARIABLE:
				result.push_back(v[i]);
				break;
			case NUMBER:
				result.push_back(v[i]);
				break;
			case BINARY_OPERATOR:
				p = priority(v[i]);
				if (v[i].s == "^") {
					while (!st.empty() && priority(st.top()) > p) {
						result.push_back(st.top());
						st.pop();
					}
					st.push(v[i]);
				} else {
					while (!st.empty() && priority(st.top()) >= p) {
						result.push_back(st.top());
						st.pop();
					}
					st.push(v[i]);
				}
				break;
			case UNARY_OPERATOR:
				p = priority(v[i]);
				while (!st.empty() && priority(st.top()) > p) {
					result.push_back(st.top());
					st.pop();
				}
				st.push(v[i]);				
				break;
			case FUNCTION:
				p = priority(v[i]);
				while (!st.empty() && priority(st.top()) > p) {
					result.push_back(st.top());
					st.pop();
				}
				st.push(v[i]);
				break;
			default:
				throw ExpressionWrongFormat();
		}
	}
	while (!st.empty()) {
		result.push_back(st.top());
		st.pop();
	}
}

// Calculate the expression result
double Expression::calculate(const vector < Token > &v) const {
	stack < double > st;
	double cur, cur2, tmp;
	map < string, double > :: const_iterator it;
	for (size_t i = 0;i < v.size();i++) {
		switch (v[i].type) {
			case NUMBER:
				sscanf(v[i].s.c_str(), "%50lg", &cur);
				st.push(cur);
				break;
			case VARIABLE:
				it = variables.find(v[i].s);
				if (it == variables.end()) {
					throw ExpressionUnknownVariable(v[i].s);
				}
				st.push(it->second);
				break;
			case UNARY_OPERATOR:
				if (st.empty()) {
					throw ExpressionWrongFormat();
				}
				cur = st.top();
				st.pop();
				st.push(applyUnary(v[i].s, cur));
				break;
			case BINARY_OPERATOR:
				if (st.size() < 2) {
					throw ExpressionWrongFormat();
				}
				cur = st.top();
				st.pop();
				cur2 = st.top();
				st.pop();
				tmp = applyBinary(v[i].s, cur2, cur);
				st.push(tmp);
				break;
			case FUNCTION:
				if (st.empty()) {
					throw ExpressionWrongFormat();
				}
				cur = st.top();
				st.pop();
				tmp = applyFunction(v[i].s, cur);
				st.push(tmp);
				break;
			default:
				throw ExpressionWrongFormat();
		}
	}
	if (st.size() != 1) {
		throw ExpressionWrongFormat();
	}
	double result = st.top();
	if (result != result || fabs(result) == numeric_limits<double>::infinity()) {
		throw ExpressionBadResult();
	}
	return result;
}

// Check if expression is valid
void Expression::checkExpression(const vector < Token > &v) const {
	int n = static_cast<int>(v.size());
	bool fail = false;

	char possibleNeighbours[][8] =  {{0, 0, 0, 0, 0, 0, 0, 0},
	                                 {0, 1, 1, 1, 1, 1, 1, 1},
	                                 {0, 1, 0, 0, 0, 1, 1, 1},
	                                 {0, 1, 0, 0, 0, 1, 1, 1},
	                                 {0, 1, 0, 0, 0, 1, 1, 1},
	                                 {0, 1, 1, 1, 0, 0, 0, 0},
	                                 {0, 1, 1, 1, 0, 0, 0, 0},
	                                 {0, 1, 0, 0, 0, 0, 0, 0}};

	fail |= (n == 0);
	for (int i = 0;i < n - 1;i++) {
		fail |= !possibleNeighbours[v[i].type][v[i + 1].type];
		fail |= (v[i].type == FUNCTION && v[i + 1].s != "(");
		fail |= (v[i].s == "(" && v[i + 1].s == ")");
		fail |= (v[i].type == UNDEFINED_OPERATOR && v[i + 1].s == ")");
		if (v[i].s == ")") {
			fail |= (v[i + 1].type != UNDEFINED_OPERATOR && v[i + 1].s != ")");
		}
		if (v[i + 1].s == "(") {
			fail |= (v[i].type != FUNCTION && v[i].type != UNDEFINED_OPERATOR && v[i].s != "(");
		}
	}
	fail |= (n >= 1 && v[n - 1].type == FUNCTION);
	fail |= (n >= 1 && v[n - 1].type == UNDEFINED_OPERATOR);

	int parenthesisBalance = 0;
	for (int i = 0;i < n;i++) {
		if (v[i].s == "(") {
			parenthesisBalance++;
		} else if (v[i].s == ")") {
			parenthesisBalance--;
		}
		fail |= (parenthesisBalance < 0);
	}
	fail |= (parenthesisBalance != 0);

	if (fail) {
		throw ExpressionWrongFormat();
	}
}

Expression::Expression(const string& s) {
	vector < Token > tokens;
	split(s, tokens);
	checkExpression(tokens);
	convertExpression(tokens, expr);
}

double Expression::eval(const string& vars) {
	setVariables(vars);
	return calculate(expr);
}
