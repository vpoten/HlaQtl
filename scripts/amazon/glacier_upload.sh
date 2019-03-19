#!/bin/bash

jardir=/path/to/lib

jars=$jardir/aws-java-sdk-1.3.33.jar:$jardir/aspectjrt.jar:$jardir/aspectjweaver.jar:$jardir/commons-codec-1.3.jar
jars=$jars:$jardir/commons-logging-1.1.1.jar:$jardir/freemarker-2.3.18.jar:$jardir/httpclient-4.1.1.jar:$jardir/httpcore-4.1.jar
jars=$jars:$jardir/jackson-core-asl-1.8.7.jar:$jardir/jackson-mapper-asl-1.8.7.jar:$jardir/mail-1.4.3.jar:$jardir/spring-beans-3.0.7.jar
jars=$jars:$jardir/spring-context-3.0.7.jar:$jardir/spring-core-3.0.7.jar:$jardir/stax-api-1.0.1.jar:$jardir/stax-1.2.0.jar

script=glacier_upload.groovy

echo groovy -cp $jars $script $1 $2
