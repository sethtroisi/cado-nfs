BEGIN{
FS=":"
}
{
  if (substr($0, 1, 1) == "#")
    {
      printf("%s\n", $0);
    }
  else
    {
      printf("%s:%s:%s\n", $1, $3, $2);
    }
}
END{
}
