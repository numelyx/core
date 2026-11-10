#include <numelyx/matrix.h>

int main(void)
{
    numelyx_matrix *identity = NULL;
    numelyx_matrix *view = NULL;
    numelyx_storage_kind kind = NUMELYX_STORAGE_DENSE;
    numelyx_matrix_flags flags = {0};

    if (numelyx_matrix_create_identity(3, &identity) != NUMELYX_SUCCESS)
        return 1;
    if (numelyx_matrix_get_storage_kind(identity, &kind) != NUMELYX_SUCCESS
        || kind != NUMELYX_STORAGE_IDENTITY)
        goto fail;
    if (numelyx_matrix_get_flags(identity, &flags) != NUMELYX_SUCCESS
        || !flags.square || !flags.diagonal || !flags.upper_triangular
        || !flags.lower_triangular || !flags.symmetric || !flags.identity
        || flags.zero)
        goto fail;
    if (numelyx_matrix_transpose(identity, &view) != NUMELYX_SUCCESS)
        goto fail;
    if (numelyx_matrix_get_storage_kind(view, &kind) != NUMELYX_SUCCESS
        || kind != NUMELYX_STORAGE_IDENTITY)
        goto fail;

    numelyx_matrix_release(view);
    numelyx_matrix_release(identity);
    return 0;

fail:
    numelyx_matrix_release(view);
    numelyx_matrix_release(identity);
    return 2;
}
