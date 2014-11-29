#include "cado.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "select_mpi.h"
#include "params.h"
#include "bw-common.h"
#include "filenames.h"
#include "portability.h"
#include "misc.h"

int command_argc;
char ** command_argv;

void * shell_prog(parallelizing_info_ptr pi, param_list pl MAYBE_UNUSED, void * arg MAYBE_UNUSED)
{
    char ** argv = malloc((command_argc + 5) * sizeof(char *));
    int i;
    size_t len = 1;
    for(i = 0 ; i < command_argc ; i++) {
        argv[i] = strdup(command_argv[i]);
        len += 1 + strlen(argv[i]);
    }
    int rc;
    rc = asprintf(&(argv[i]), "%d", pi->wr[0]->jrank);
    ASSERT_ALWAYS(rc != -1);
    len += 1 + strlen(argv[i]);
    i++;
    rc = asprintf(&(argv[i]), "%d", pi->wr[1]->jrank);
    ASSERT_ALWAYS(rc != -1);
    len += 1 + strlen(argv[i]);
    i++;
    rc = asprintf(&(argv[i]), "%d", pi->wr[0]->trank);
    ASSERT_ALWAYS(rc != -1);
    len += 1 + strlen(argv[i]);
    i++;
    rc = asprintf(&(argv[i]), "%d", pi->wr[1]->trank);
    ASSERT_ALWAYS(rc != -1);
    len += 1 + strlen(argv[i]);
    i++;
    argv[i++]=NULL;

    char * cmdline = malloc(len);
    size_t k = 0;
    for(i = 0 ; i < command_argc + 4 ; i++) {
        snprintf(cmdline + k, len-k, " %s", argv[i]);
        k += 1 + strlen(argv[i]);
    }

    pid_t child = fork();


    if (child == 0) {
        const char * script = param_list_lookup_string(pl, "script");
        const char * command = param_list_lookup_string(pl, "command");
        if (script) {
            execvp(script, argv);
        } else if (command) {
            int rc = system(command);
            exit(WEXITSTATUS(rc));
        } else {
            fprintf(stderr, "Please supply either script= or command=\n");
            exit(1);
        }
    } else {
        int status = 0;
        waitpid(child,&status,0);
        if (WIFEXITED(status)) {
            int rc;
            if ((rc = WEXITSTATUS(status)) != 0) {
                fprintf(stderr, "Command%s exited with status %d\n", cmdline, rc);
            }
        } else if (WIFSIGNALED(status)) {
            fprintf(stderr, "Command%s exited with signal %d\n", cmdline, WTERMSIG(status));
        } else {
            fprintf(stderr, "Command%s fooleed wait() !\n", cmdline);
        }
    }
    serialize(pi->m);
    for(i = 0 ; i < command_argc + 4 ; i++) {
        free(argv[i]);
    }
    free(argv);
    free(cmdline);
    return NULL;
}

int main(int argc, char * argv[])
{
    param_list pl;
    bw_common_init_new(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    param_list_add_key(pl, "mn", "0", PARAMETER_FROM_FILE);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);
    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);

    param_list_warn_unused(pl);
    
    /* not sure what I should do with this
    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    } */

    command_argc = argc;
    command_argv = argv;

    pi_go(shell_prog, pl, 0);

    param_list_clear(pl);
    bw_common_clear_new(bw);

    return 0;
}

