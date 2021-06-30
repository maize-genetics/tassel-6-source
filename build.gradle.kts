import org.gradle.api.JavaVersion.VERSION_11
import org.jetbrains.dokka.gradle.DokkaTask
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

// Note Kotlin version needs to be updated in both the buildscript and plugins.
// Dependencies will follow the buildscript

group = "net.maizegenetics"
version = "6.0"

buildscript {
    val kotlinVersion by extra ("1.4.32")

    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven("https://dl.bintray.com/kotlin/kotlin-eap")
        maven("https://plugins.gradle.org/m2/")
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
        classpath("org.jetbrains.dokka:dokka-gradle-plugin:1.4.30")
    }
}


plugins {

    val kotlinVersion = "1.4.32"
    java
    kotlin("jvm") version kotlinVersion
    kotlin("plugin.serialization") version kotlinVersion

    // Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "5.2.0"
    id("application")
    id("org.openjfx.javafxplugin") version "0.0.10"

    id("org.jetbrains.dokka") version "1.4.30"
    `java-library`
    `maven-publish`
    signing

}

application {
    mainClassName = "net.maizegenetics.pipeline.TasselPipeline"
}

apply {
    plugin("kotlinx-serialization")
    plugin("org.jetbrains.dokka")
}

repositories {
    mavenCentral()
    jcenter()
    maven("https://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
    maven("https://dl.bintray.com/kotlin/kotlin-eap")
    maven("https://kotlin.bintray.com/kotlinx")
    maven("https://oss.sonatype.org/content/repositories/snapshots/")
}

javafx {
    modules("javafx.controls", "javafx.web")
}

dependencies {

    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("log4j:log4j:1.2.17")
    implementation("com.google.guava:guava:29.0-jre")
    implementation("com.github.samtools:htsjdk:2.24.1")
    implementation("com.googlecode.efficient-java-matrix-library:ejml:0.23")
    implementation("colt:colt:1.2.0")
    implementation("com.googlecode.json-simple:json-simple:1.1.1")
    implementation("org.glassfish:javax.json:1.1.4")
    implementation("commons-codec:commons-codec:1.11")
    implementation("org.apache.commons:commons-math3:3.6.1")

    implementation("org.openjfx:javafx-controls:11.0.2")
    implementation("org.openjfx:javafx-web:11.0.2")

    implementation("org.jetbrains.kotlin:kotlin-stdlib:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-stdlib-common:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-stdlib-jdk8:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-stdlib-jdk7:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.4.3")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-core:1.1.0")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.1.0")

    implementation("khttp:khttp:1.0.0")

    val kotestVersion = "4.2.6"
    listOf("runner-junit5", "assertions-core", "property").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }

}

java {
    sourceCompatibility = VERSION_11
    targetCompatibility = VERSION_11
    withSourcesJar()
}

tasks.withType<KotlinCompile>().configureEach {
    kotlinOptions.jvmTarget = "11"
}

tasks {
    println("Source directories: ${sourceSets["main"].allSource.srcDirs}")
}

tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
}

val dokkaHtml by tasks.getting(org.jetbrains.dokka.gradle.DokkaTask::class)

val dokkaJar by tasks.creating(Jar::class) {
    dependsOn(dokkaHtml)
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "TASSEL 6: ${property("version")}"
    archiveClassifier.set("javadoc")
    from(dokkaHtml.outputDirectory)
}

publishing {
    publications {

        create<MavenPublication>("maven") {
            artifactId = "tassel6"
            description = "net.maizegenetics:tassel6:$version"

            from(components["java"])
            artifact(dokkaJar)

            versionMapping {
                usage("java-api") {
                    fromResolutionOf("runtimeClasspath")
                }
                usage("java-runtime") {
                    fromResolutionResult()
                }
            }

            repositories {
                maven {
                    val releasesRepoUrl = "https://oss.sonatype.org/service/local/staging/deploy/maven2"
                    val snapshotsRepoUrl = "https://oss.sonatype.org/content/repositories/snapshots"
                    url = uri(if (version.toString().endsWith("SNAPSHOT")) snapshotsRepoUrl else releasesRepoUrl)
                    credentials {
                        try {
                            username = property("ossrhUsername") as String?
                            password = property("ossrhPassword") as String?
                        } catch (e: Exception) {
                            println("Unable to get username and password for nexus maven central release.")
                        }
                    }
                }
            }

            pom {
                name.set("tassel-6")
                description.set("TASSEL 6 is a software package to evaluate traits association. Feature Tables are at the heart of the package where, a feature is a range of positions or a single position. Row in the that table are taxon.")
                url.set("http://www.maizegenetics.net/tassel")
                licenses {
                    license {
                        name.set("The Apache License, Version 2.0")
                        url.set("http://www.apache.org/licenses/LICENSE-2.0.txt")
                    }
                }
                developers {
                    developer {
                        name.set("Ed Buckler")
                        email.set("esb33@cornell.edu")
                    }
                    developer {
                        id.set("tmc46")
                        name.set("Terry Casstevens")
                        email.set("tmc46@cornell.edu")
                    }
                    developer {
                        name.set("Zack Miller")
                        email.set("zrm22@cornell.edu")
                    }
                    developer {
                        name.set("Lynn Johnson")
                        email.set("lcj34@cornell.edu")
                    }
                    developer {
                        name.set("Brandon Monier")
                        email.set("bm646@cornell.edu")
                    }
                    developer {
                        name.set("Peter Bradbury")
                        email.set("pjb39@cornell.edu")
                    }
                }
                scm {
                    connection.set("scm:git:git://bitbucket.org:tasseladmin/tassel-6-source.git")
                    developerConnection.set("scm:git:ssh://bitbucket.org:tasseladmin/tassel-6-source.git</developer")
                    url.set("https://bitbucket.org/tasseladmin/tassel-6-source/src")
                }
            }
        }
    }
}

signing {
    useGpgCmd()
    sign(publishing.publications["maven"])
}

tasks.javadoc {
    dependsOn("dokkaJavadoc")
    if (JavaVersion.current().isJava9Compatible) {
        (options as StandardJavadocDocletOptions).addBooleanOption("html5", true)
    }
}

