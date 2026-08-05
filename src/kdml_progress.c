#include "kdml.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _WIN32
# include <windows.h>
#else
# include <sys/time.h>
#endif

#define KDML_PROGRESS_BAR_WIDTH 32
#define KDML_PROGRESS_CHECKPOINTS 200
#define KDML_PROGRESS_MIN_SECONDS 0.5

static double kdml_progress_now(void)
{
#ifdef _WIN32
    LARGE_INTEGER counter;
    LARGE_INTEGER frequency;

    if (QueryPerformanceFrequency(&frequency) &&
        QueryPerformanceCounter(&counter) &&
        frequency.QuadPart > 0) {
        return (double) counter.QuadPart / (double) frequency.QuadPart;
    }
    return (double) GetTickCount64() / 1000.0;
#else
    {
        struct timeval value;
        if (gettimeofday(&value, NULL) == 0) {
            return (double) value.tv_sec + (double) value.tv_usec / 1e6;
        }
    }
    return 0.0;
#endif
}

static void kdml_progress_duration(double seconds, char *buffer,
                                   size_t buffer_size)
{
    unsigned long long whole;
    unsigned long long hours;
    unsigned long long minutes;

    if (!R_FINITE(seconds) || seconds < 0.0) {
        seconds = 0.0;
    }
    whole = (unsigned long long) floor(seconds);
    hours = whole / 3600ULL;
    minutes = (whole % 3600ULL) / 60ULL;
    whole %= 60ULL;
    snprintf(buffer, buffer_size, "%02llu:%02llu:%02llu",
             hours, minutes, whole);
}

static void kdml_progress_format(const kdml_progress *progress,
                                 int chain, long long completed,
                                 const char *status, double now,
                                 char *line, size_t line_size)
{
    char bar[KDML_PROGRESS_BAR_WIDTH + 1];
    char elapsed_text[32];
    char eta_text[32];
    char rate_text[64];
    char label[64];
    double elapsed;
    double fraction;
    double eta;
    int filled;

    if (now < progress->start_seconds) {
        now = progress->start_seconds;
    }
    elapsed = now - progress->start_seconds;
    fraction = progress->total > 0 ?
        (double) completed / (double) progress->total : 1.0;
    if (fraction < 0.0) {
        fraction = 0.0;
    } else if (fraction > 1.0) {
        fraction = 1.0;
    }
    filled = (int) floor(
        fraction * (double) KDML_PROGRESS_BAR_WIDTH + 0.5
    );
    memset(bar, '#', (size_t) filled);
    memset(bar + filled, '-',
           (size_t) (KDML_PROGRESS_BAR_WIDTH - filled));
    bar[KDML_PROGRESS_BAR_WIDTH] = '\0';

    kdml_progress_duration(elapsed, elapsed_text, sizeof(elapsed_text));
    if (completed <= 0) {
        snprintf(eta_text, sizeof(eta_text), "--:--:--");
        snprintf(rate_text, sizeof(rate_text), "--.--it/s");
    } else {
        if (completed >= progress->total) {
            eta = 0.0;
        } else {
            eta = elapsed *
                (double) (progress->total - completed) /
                (double) completed;
        }
        kdml_progress_duration(eta, eta_text, sizeof(eta_text));
        if (elapsed > 0.0) {
            snprintf(rate_text, sizeof(rate_text), "%.2fit/s",
                     (double) completed / elapsed);
        } else {
            snprintf(rate_text, sizeof(rate_text), "--.--it/s");
        }
    }

    if (chain == 0) {
        snprintf(label, sizeof(label), "main chains");
    } else if (progress->chains == 1) {
        snprintf(label, sizeof(label), "main chain");
    } else {
        snprintf(label, sizeof(label), "main chain %d/%d",
                 chain, progress->chains);
    }
    snprintf(
        line, line_size,
        "%s: %3.0f%%|%s| %lld/%lld [%s<%s, %7s] %s",
        label, 100.0 * fraction, bar,
        completed, progress->total, elapsed_text, eta_text,
        rate_text, status
    );
}

static void kdml_progress_render(kdml_progress *progress,
                                 const char *status, double now)
{
    char line[256];
    char padding[256];

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    kdml_progress_format(
        progress, progress->aggregate ? 0 : progress->chain,
        progress->completed, status, now, line, sizeof(line)
    );
    if (progress->in_place) {
        const int width = (int) strlen(line);
        int extra = progress->last_width - width;
        if (extra < 0) {
            extra = 0;
        } else if (extra >= (int) sizeof(padding)) {
            extra = (int) sizeof(padding) - 1;
        }
        memset(padding, ' ', (size_t) extra);
        padding[extra] = '\0';
        REprintf("\r%s%s", line, padding);
        progress->last_width = width;
    } else {
        REprintf("%s\n", line);
    }
    R_FlushConsole();
}

void kdml_progress_initialize(kdml_progress *progress, int enabled,
                              int in_place, int chains, long long total)
{
    if (progress == NULL) {
        return;
    }
    memset(progress, 0, sizeof(*progress));
    progress->enabled = enabled != 0;
    progress->in_place = in_place != 0;
    progress->chains = chains;
    progress->total = total;
    progress->stride = total / KDML_PROGRESS_CHECKPOINTS;
    if (progress->stride < 1) {
        progress->stride = 1;
    }
}

void kdml_progress_begin(kdml_progress *progress, int chain)
{
    double now;

    if (progress == NULL || !progress->enabled) {
        return;
    }
    now = kdml_progress_now();
    progress->active = 1;
    progress->chain = chain;
    progress->aggregate = chain == 0;
    progress->completed = 0;
    progress->last_width = 0;
    progress->start_seconds = now;
    progress->last_emit_seconds = now;
    kdml_progress_render(progress, "running", now);
}

void kdml_progress_tick(kdml_progress *progress, long long completed)
{
    double now;

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    if (completed < 0) {
        completed = 0;
    } else if (completed > progress->total) {
        completed = progress->total;
    }
    progress->completed = completed;
    if (completed >= progress->total ||
        (completed != 1 && completed % progress->stride != 0)) {
        return;
    }
    now = kdml_progress_now();
    if (now - progress->last_emit_seconds < KDML_PROGRESS_MIN_SECONDS) {
        return;
    }
    kdml_progress_render(progress, "running", now);
    progress->last_emit_seconds = now;
}

void kdml_progress_end(kdml_progress *progress, const char *status)
{
    double now;

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    if (status == NULL) {
        status = "failed";
    }
    now = kdml_progress_now();
    kdml_progress_render(progress, status, now);
    if (progress->in_place) {
        REprintf("\n");
        R_FlushConsole();
    }
    progress->last_emit_seconds = now;
    progress->active = 0;
}

static void kdml_progress_group_render(kdml_progress *progress,
                                       const int *completed,
                                       const char *status, double now,
                                       int redraw)
{
    int chain;

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    if (progress->in_place && redraw) {
        REprintf("\033[%dA", progress->chains);
    }
    for (chain = 0; chain < progress->chains; ++chain) {
        char line[256];
        const char *chain_status = status;
        long long value = completed == NULL ? 0 : completed[chain];
        if (value < 0) {
            value = 0;
        } else if (value > progress->total) {
            value = progress->total;
        }
        if (value >= progress->total &&
            (strcmp(status, "running") == 0 ||
             strcmp(status, "failed") == 0 ||
             strcmp(status, "interrupted") == 0)) {
            chain_status = "done";
        }
        kdml_progress_format(
            progress, chain + 1, value, chain_status, now,
            line, sizeof(line)
        );
        if (progress->in_place) {
            REprintf("\r%s\033[K\n", line);
        } else {
            REprintf("%s\n", line);
        }
    }
    R_FlushConsole();
}

void kdml_progress_group_begin(kdml_progress *progress,
                               const int *completed)
{
    double now;

    if (progress == NULL || !progress->enabled) {
        return;
    }
    now = kdml_progress_now();
    progress->active = 1;
    progress->aggregate = 0;
    progress->start_seconds = now;
    progress->last_emit_seconds = now;
    kdml_progress_group_render(
        progress, completed, "running", now, 0
    );
}

void kdml_progress_group_tick(kdml_progress *progress,
                              const int *completed)
{
    double now;

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    now = kdml_progress_now();
    if (now - progress->last_emit_seconds < KDML_PROGRESS_MIN_SECONDS) {
        return;
    }
    kdml_progress_group_render(
        progress, completed, "running", now, 1
    );
    progress->last_emit_seconds = now;
}

void kdml_progress_group_end(kdml_progress *progress,
                             const int *completed, const char *status)
{
    double now;

    if (progress == NULL || !progress->enabled || !progress->active) {
        return;
    }
    if (status == NULL) {
        status = "failed";
    }
    now = kdml_progress_now();
    kdml_progress_group_render(
        progress, completed, status, now, 1
    );
    progress->last_emit_seconds = now;
    progress->active = 0;
}
