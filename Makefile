# Makefile for the 'STORM' program

NAME = RegExpMatch
LIBS =  -lpopt -lm  -lboost_regex#-lefence
CC = gcc
CPP = g++
DEBUGFLAGS = -Wall -g
CFLAGS = #-static #-O2
MAIN = main.cpp

$(NAME):	$(MAIN)
	$(CPP) $(DEBUGFLAGS) $(CFLAGS) -o $(NAME) $(MAIN) $(LIBS)

