// Math expression parser
// @author Stepan Belousov
#pragma once
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

// Main class for a math expression
class Expression {
	static const double pi = M_PI;

	enum TokenType { UNDEFINED, PARENTHESIS, UNDEFINED_OPERATOR, 
	BINARY_OPERATOR, UNARY_OPERATOR, VARIABLE, NUMBER, FUNCTION };

	static const std::string TokenTypeList[];

	std::set < std::string > functions;
	std::set < std::string > operators;
	std::map < std::string, double > variables;

	// Service class for one token
	class Token {
	public:
		std::string s;
		TokenType type;
		Token(const std::string& s, const TokenType& type = UNDEFINED): s(s), type(type) {}
		void print() const {
			printf("%s   <%s>\n", s.c_str(), TokenTypeList[type].c_str());
		}
	};

	std::vector < Token > expr;

	// Implementation functions
	void initialize();
	double applyFunction(const std::string& s, double arg) const;
	double applyUnary(const std::string& s, double arg) const;
	double applyBinary(const std::string& s, double arg1, double arg2) const;
	bool isOperator(const std::string& s) const;
	bool isDelimiter(char c) const;
	bool isVariable(const std::string& s) const;
	TokenType getType(const Token& t) const;
	void markUnaries(std::vector < Token > & v) const;
	void split(const std::string& s, std::vector < Token > & result);
	int priority(const Token& t) const;
	void convertExpression(std::vector < Token > & v, std::vector < Token > & result) const;
	double calculate(const std::vector < Token > & v) const;
	void checkExpression(const std::vector < Token > & v) const;
	void setVariables(const std::string& s);

// Public interface
public:
	Expression(const std::string& s);
	double eval(const std::string& vars = "");
};

class ExpressionWrongFormat : public std::exception {
public:
	const char* what() const throw() {
		return "Input expression is not in correct format";
	}
};

class ExpressionBadResult : public std::exception {
public:
	const char* what() const throw() {
		return "Expression result is either undefined or infinite";
	}
};

class ExpressionUnknownVariable : public std::exception {
	std::string message;
public:
	ExpressionUnknownVariable(const std::string& name): message("Unknown variable: " + name) {}
	const char* what() const throw() {
		return message.c_str();
	}
	~ExpressionUnknownVariable() throw() {}
};
