if status is-interactive
    # Commands to run in interactive sessions can go here
end

function fish_greeting; end

alias ls='exa -la'
set -x TERM xterm-256color
