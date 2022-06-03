_cg_spica () {
  local args
  local c=${COMP_WORDS[COMP_CWORD]}
  local p=${COMP_WORDS[1]}
  if [ $COMP_CWORD -le 1 ] ;then
      args="json2top map2cg maptraj ENM"
      COMPREPLY=( `compgen -W "$args" -- $c`)
      return 0
  fi
  case $p in 
  json2top )
      :;;
  map2cg )
      COMPREPLY=( `compgen -S ' ' -f -- $c; compgen -S '/' -d $c`);;
  maptraj )
      COMPREPLY=( `compgen -S ' ' -f -- $c; compgen -S '/' -d $c`);;
  ENM )
      COMPREPLY=( `compgen -S ' ' -f -- $c; compgen -S '/' -d $c`);;
  esac
  #COMPREPLY=( `compgen -W "$args" -- $c; compgen -S '/' -d $c`)
}
complete -F _cg_spica cg_spica
