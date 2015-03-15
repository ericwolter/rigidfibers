#ifndef FIBERS_FIBEROPT_H_
#define FIBERS_FIBEROPT_H_
/*
 *  fiberopt.h - command line arguments parsing
 *
 *  Copyright (C) 2014  Eric Wolter <eric.wolter@gmx.de>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>

typedef struct
{
    /* special */
    const char *usage_pattern;
    const char *help_message;
    /* commands */

    /* arguments */
    char *layout;
    /* options without arguments */
    int gui;
    int help;
    int version;
    /* options with arguments */
} FiberArgs;

const char help_message[] =
    "Rigid Fibers.\n"
    "\n"
    "Usage:\n"
    "  fibers <layout>\n"
    "  fibers --help\n"
    "  fibers --version\n"
    "\n"
    "Options:\n"
    "  -h --help     Show this screen.\n"
    "  -v --version  Show version.\n"
    "\n"
    "";

const char usage_pattern[] =
    "Usage:\n"
    "Usage:\n"
    "  fibers <layout>\n"
    "  fibers --help\n"
    "  fibers --version";

typedef struct
{
    const char *name;
    bool value;
    char pad[7];            // nothing we can do about this padding so at least
    // make it explicit
} Command;

typedef struct
{
    const char *name;       // the name of the argument, only used for help
    char *value;            // the value of the argument
} Argument;

typedef struct
{
    const char *oshort;     // a short name for this option, e.g. -h
    const char *olong;      // a long name for this option, e.g. --help
    char *argument;         // a potential list of arguments for this option
    bool argcount;          // how many arguments does this option have
    bool value;             // a flag indicating that this option was supplied
    char pad[6];            // nothing we can do about this padding so at least
    // make it explicit
} Option;

// Describes the overall structure of the arguments
typedef struct
{
    Command *commands;      // pointer to all command structs
    Argument *arguments;    // pointer to all arugment structs
    Option *options;        // pointer to all options structs
    int n_commands;         // how many commands are there in total?
    int n_arguments;        // how many arguments are there in total?
    int n_options;          // how many options are there in total?
    char pad[4];            // nothing we can do about this padding so at least
    // make it explicit
} Elements;

// Input arguments are bissected into tokens to be parse one after the other
typedef struct
{
    int i;                  // the current argument index
    int argc;               // the total number of arguments
    char *current;          // the current token to be parsed
    char **argv;            // the input arguments as supplied from the
    // commandline
} Tokens;

// The prototypes declaring all the functions needed for command line parsing
Tokens tokens_new(int argc, char **argv);
Tokens *tokens_move(Tokens *ts);
int parse_args(Tokens *ts, Elements *elements);
int elems_to_args(Elements *elements, FiberArgs *args, bool help,
                  const char *version);
int parse_argcmd(Tokens *ts, Elements *elements);
int parse_long(Tokens *ts, Elements *elements);
int parse_short(Tokens *ts, Elements *elements);
FiberArgs fiberopt(int argc, char *argv[], bool help, const char *version);


// Starts the parsing of the input arguments as supplied from the commandline
Tokens tokens_new(int argc, char **argv)
{
    char* firstArg = NULL;
    if (argc > 1)
    {
        firstArg = argv[1];
    }

    // ignore the program name start with argument at index 1
    Tokens ts = {1, argc - 1, firstArg, argv};
    return ts;
}

Tokens *tokens_move(Tokens *ts)
{
    if (ts->i < ts->argc)
    {
        ts->current = ts->argv[++ts->i];
    }
    else
    {
        ts->current = NULL;
    }
    return ts;
}

int parse_argcmd(Tokens *ts, Elements *elements)
{
    // currently we only need arguments
    //int n_commands = elements->n_commands;
    int n_arguments = elements->n_arguments;
    Argument *argument;
    Argument *arguments = elements->arguments;

    if (ts->i > n_arguments)
    {
        fprintf(stderr, "too many arguments\n");
        return 1;
    }

    argument = &arguments[ts->i - 1];
    argument->value = ts->current;

    tokens_move(ts);

    return 0;
}

int parse_long(Tokens *ts, Elements *elements)
{
    int n_options = elements->n_options;

    // look for option with argument
    char *eq = strchr(ts->current, '=');

    Option *option = NULL;
    Option *options = elements->options;

    size_t len_prefix = (unsigned long)(eq - (ts->current)) / sizeof(char);

    int i;
    for (i = 0; i < n_options; i++)
    {
        option = &options[i];
        if (!strncmp(ts->current, option->olong, len_prefix))
        {
            break;
        }
    }
    if (i == n_options)
    {
        fprintf(stderr, "%s is not recognized\n", ts->current);
    }
    tokens_move(ts);
    if (option->argcount)
    {

    }
    else
    {
        if (eq != NULL)
        {
            fprintf(stderr, "%s must not have an argument\n", option->olong);
            return 1;
        }
        option->value = true;
    }
    return 0;
}

int parse_short(Tokens *ts, Elements *elements)
{
    (void)(ts);
    (void)(elements);
    return 1;
}

int parse_args(Tokens *ts, Elements *elements)
{
    int ret;

    while (ts->current != NULL)
    {
        if (ts->current[0] == '-' && ts->current[1] == '-')
        {
            ret = parse_long(ts, elements);
        }
        else if (ts->current[0] == '-' && ts->current[1] != '\0')
        {
            ret = parse_short(ts, elements);
        }
        else
        {
            ret = parse_argcmd(ts, elements);
        }
        if (ret) return ret;
    }
    return 0;
}

int elems_to_args(Elements *elements, FiberArgs *args, bool help,
                  const char *version)
{
    // Command *command;
    Argument *argument;
    Option *option;
    int i;

    /* options */
    for (i = 0; i < elements->n_options; i++)
    {
        option = &elements->options[i];
        if (help && option->value && !strcmp(option->olong, "--help"))
        {
            printf("%s", args->help_message);
            return 1;
        }
        else if (version && option->value &&
                 !strcmp(option->olong, "--version"))
        {
            printf("%s\n", version);
            return 1;
        }
        else if (!strcmp(option->olong, "--gui"))
        {
            args->gui = option->value;
        }
    }
    // /* commands */
    // for (i=0; i < elements->n_commands; i++) {
    //     command = &elements->commands[i];
    //     if (!strcmp(command->name, "create")) {
    //         args->create = command->value;
    //     }
    // }
    /* arguments */
    for (i = 0; i < elements->n_arguments; i++)
    {
        argument = &elements->arguments[i];
        if (!strcmp(argument->name, "<layout>"))
        {
            args->layout = argument->value;
        }
    }

    // both arguments are required
    if (!args->layout)
    {
        fprintf(stderr, "argument <layout> is required\n");
    }
    if (!args->layout)
    {
        exit(EXIT_FAILURE);
    }

    return 0;
}

FiberArgs fiberopt(int argc, char *argv[], bool help, const char *version)
{
    FiberArgs args =
    {
        usage_pattern, help_message, NULL, 0, 0, 0
    };
    Command *commands = NULL;
    Argument arguments[] =
    {
        {"<layout>", NULL}
    };
    Option options[] =
    {
        {NULL, "--gui", NULL, 0, 0, {0}},
        {"-h", "--help", NULL, 0, 0, {0}},
        {"-v", "--version", NULL, 0, 0, {0}}
    };
    Elements elements = {commands, arguments, options, 0, 1, 3, {0}};

    Tokens ts = tokens_new(argc, argv);
    if (parse_args(&ts, &elements))
    {
        exit(EXIT_FAILURE);
    }
    if (elems_to_args(&elements, &args, help, version))
    {
        exit(EXIT_SUCCESS);
    }
    return args;
}

#endif // FIBERS_FIBEROPT_H_
