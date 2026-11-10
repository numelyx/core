#ifndef NUMELYX_EXPORT_H
#define NUMELYX_EXPORT_H

#if defined(_WIN32) && defined(NUMELYX_SHARED)
#    if defined(NUMELYX_BUILDING_LIBRARY)
#        define NUMELYX_API __declspec(dllexport)
#    else
#        define NUMELYX_API __declspec(dllimport)
#    endif
#elif defined(__GNUC__) && (__GNUC__ >= 4)
#    define NUMELYX_API __attribute__((visibility("default")))
#else
#    define NUMELYX_API
#endif

#ifdef __cplusplus
#    define NUMELYX_BEGIN_DECLS extern "C" {
#    define NUMELYX_END_DECLS }
#else
#    define NUMELYX_BEGIN_DECLS
#    define NUMELYX_END_DECLS
#endif

#endif /* NUMELYX_EXPORT_H */
