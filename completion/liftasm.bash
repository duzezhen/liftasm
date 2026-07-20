_liftasm_completion()
{
    local cur
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"

    if [[ $COMP_CWORD -eq 1 ]]; then
        local subcommands="
            stat
            seq
            gfa2fa
            depth
            bubble
            augment
            clean
            collapse
            file2map
            liftover
            mapq_boost
            deoverlap
            align
            ressit
            split
        "

        COMPREPLY=($(compgen -W "$subcommands" -- "$cur"))
        return 0
    fi

    COMPREPLY=($(compgen -f -- "$cur"))
    return 0
}


_liftasm_register_completion()
{
    local cmd

    while IFS= read -r cmd; do
        [[ "$cmd" == liftasm* ]] || continue
        complete -F _liftasm_completion "$cmd"
    done < <(compgen -c | sort -u)
}

_liftasm_register_completion
