package datagenerators;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

/**
 * Turns a packet trace into the stream format of the experiments (segments mode): one packet record
 * per line, {@code src,dst,proto}, in capture order. This is the symbolization of the thesis
 * (Sec. 5.1, the 3-tuple ⟨ip1, ip2, proto⟩); the index hashes each line into one symbol.
 *
 * <p>Reads classic pcap, plain or gzipped, with Ethernet (802.1Q/802.1ad tags allowed), raw IP or
 * Linux cooked framing. IPv4 addresses are dotted quads; IPv6 addresses are eight lowercase hex
 * groups without leading zeros and without "::" compression, so every address has one spelling.
 * The protocol is the IPv4 protocol field, or for IPv6 the next header after any hop-by-hop,
 * routing, destination-options, fragment and AH headers (the last one visible if the capture is
 * truncated). Packets that are not IP, or too short to hold the addresses, are skipped and counted.
 *
 * <p>Usage: {@code PcapTokenizer <trace.pcap[.gz]> <out.txt>}
 */
public final class PcapTokenizer {

    private static final int LINK_ETHERNET = 1;
    private static final int LINK_RAW = 101;
    private static final int LINK_RAW_ALT = 12;
    private static final int LINK_LINUX_SLL = 113;

    private PcapTokenizer() {
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 2) {
            System.err.println("usage: PcapTokenizer <trace.pcap[.gz]> <out.txt>");
            System.exit(2);
        }
        Stats stats = tokenize(Path.of(args[0]), Path.of(args[1]));
        System.out.println(stats);
    }

    static Stats tokenize(Path trace, Path out) throws IOException {
        Stats stats = new Stats();
        long t0 = System.nanoTime();
        try (InputStream raw = Files.newInputStream(trace);
             InputStream unzipped = trace.toString().endsWith(".gz")
                     ? new GZIPInputStream(raw, 1 << 20) : raw;
             DataInputStream in = new DataInputStream(new BufferedInputStream(unzipped, 1 << 20));
             OutputStream os = new BufferedOutputStream(Files.newOutputStream(out), 1 << 20)) {

            byte[] global = new byte[24];
            in.readFully(global);
            int magic = readInt(global, 0, false);
            boolean swapped;
            if (magic == 0xa1b2c3d4 || magic == 0xa1b23c4d) {
                swapped = false;
            } else if (magic == 0xd4c3b2a1 || magic == 0x4d3cb2a1) {
                swapped = true;
            } else {
                throw new IOException("not a classic pcap file (magic " + Integer.toHexString(magic)
                        + "); convert pcapng with editcap -F pcap first");
            }
            int linkType = readInt(global, 20, swapped) & 0x0fffffff;
            if (linkType != LINK_ETHERNET && linkType != LINK_RAW && linkType != LINK_RAW_ALT
                    && linkType != LINK_LINUX_SLL) {
                throw new IOException("unsupported link type " + linkType);
            }

            byte[] header = new byte[16];
            byte[] packet = new byte[1 << 16];
            byte[] line = new byte[128];
            while (true) {
                try {
                    in.readFully(header);
                } catch (EOFException end) {
                    break;
                }
                int captured = readInt(header, 8, swapped);
                if (captured < 0 || captured > packet.length) {
                    throw new IOException("corrupt record after " + stats.packets + " packets (length " + captured + ")");
                }
                in.readFully(packet, 0, captured);
                stats.packets++;
                int length = format(packet, captured, linkType, line, stats);
                if (length > 0) {
                    os.write(line, 0, length);
                    stats.written++;
                }
            }
        }
        stats.seconds = (System.nanoTime() - t0) / 1e9;
        return stats;
    }

    /** Writes "src,dst,proto\n" into {@code line} and returns its length, or 0 to skip the packet. */
    private static int format(byte[] p, int len, int linkType, byte[] line, Stats stats) {
        int off;
        int etherType;
        if (linkType == LINK_ETHERNET) {
            if (len < 14) {
                stats.truncated++;
                return 0;
            }
            off = 12;
            etherType = u16(p, off);
            off += 2;
            while ((etherType == 0x8100 || etherType == 0x88a8 || etherType == 0x9100) && off + 4 <= len) {
                etherType = u16(p, off + 2);
                off += 4;
            }
        } else if (linkType == LINK_LINUX_SLL) {
            if (len < 16) {
                stats.truncated++;
                return 0;
            }
            etherType = u16(p, 14);
            off = 16;
        } else {
            if (len < 1) {
                stats.truncated++;
                return 0;
            }
            off = 0;
            int version = (p[0] & 0xff) >>> 4;
            etherType = (version == 4) ? 0x0800 : (version == 6) ? 0x86dd : -1;
        }

        int pos = 0;
        if (etherType == 0x0800) {
            if (off + 20 > len || ((p[off] & 0xff) >>> 4) != 4) {
                stats.truncated++;
                return 0;
            }
            pos = ipv4(p, off + 12, line, pos);
            line[pos++] = ',';
            pos = ipv4(p, off + 16, line, pos);
            line[pos++] = ',';
            pos = decimal(p[off + 9] & 0xff, line, pos);
            stats.ipv4++;
        } else if (etherType == 0x86dd) {
            if (off + 40 > len || ((p[off] & 0xff) >>> 4) != 6) {
                stats.truncated++;
                return 0;
            }
            pos = ipv6(p, off + 8, line, pos);
            line[pos++] = ',';
            pos = ipv6(p, off + 24, line, pos);
            line[pos++] = ',';
            pos = decimal(ipv6Protocol(p, off, len), line, pos);
            stats.ipv6++;
        } else {
            stats.nonIp++;
            return 0;
        }
        line[pos++] = '\n';
        return pos;
    }

    /** Upper-layer protocol of an IPv6 packet: skips extension headers while they are captured. */
    private static int ipv6Protocol(byte[] p, int off, int len) {
        int next = p[off + 6] & 0xff;
        int pos = off + 40;
        while (true) {
            switch (next) {
                case 0, 43, 60 -> {            // hop-by-hop, routing, destination options
                    if (pos + 2 > len) {
                        return next;
                    }
                    int skip = ((p[pos + 1] & 0xff) + 1) * 8;
                    next = p[pos] & 0xff;
                    pos += skip;
                }
                case 44 -> {                   // fragment
                    if (pos + 1 > len) {
                        return next;
                    }
                    next = p[pos] & 0xff;
                    pos += 8;
                }
                case 51 -> {                   // authentication header
                    if (pos + 2 > len) {
                        return next;
                    }
                    int skip = ((p[pos + 1] & 0xff) + 2) * 4;
                    next = p[pos] & 0xff;
                    pos += skip;
                }
                default -> {
                    return next;
                }
            }
        }
    }

    private static int ipv4(byte[] p, int at, byte[] out, int pos) {
        for (int i = 0; i < 4; i++) {
            if (i > 0) {
                out[pos++] = '.';
            }
            pos = decimal(p[at + i] & 0xff, out, pos);
        }
        return pos;
    }

    private static final byte[] HEX = "0123456789abcdef".getBytes();

    private static int ipv6(byte[] p, int at, byte[] out, int pos) {
        for (int g = 0; g < 8; g++) {
            if (g > 0) {
                out[pos++] = ':';
            }
            int v = u16(p, at + 2 * g);
            boolean started = false;
            for (int shift = 12; shift >= 0; shift -= 4) {
                int nibble = (v >>> shift) & 0xf;
                if (nibble != 0 || started || shift == 0) {
                    out[pos++] = HEX[nibble];
                    started = true;
                }
            }
        }
        return pos;
    }

    private static int decimal(int v, byte[] out, int pos) {
        if (v >= 100) {
            out[pos++] = (byte) ('0' + v / 100);
            out[pos++] = (byte) ('0' + (v / 10) % 10);
        } else if (v >= 10) {
            out[pos++] = (byte) ('0' + v / 10);
        }
        out[pos++] = (byte) ('0' + v % 10);
        return pos;
    }

    private static int u16(byte[] p, int at) {
        return ((p[at] & 0xff) << 8) | (p[at + 1] & 0xff);
    }

    private static int readInt(byte[] b, int at, boolean littleEndian) {
        if (littleEndian) {
            return (b[at] & 0xff) | (b[at + 1] & 0xff) << 8 | (b[at + 2] & 0xff) << 16 | (b[at + 3] & 0xff) << 24;
        }
        return (b[at] & 0xff) << 24 | (b[at + 1] & 0xff) << 16 | (b[at + 2] & 0xff) << 8 | (b[at + 3] & 0xff);
    }

    static final class Stats {
        long packets;
        long written;
        long ipv4;
        long ipv6;
        long nonIp;
        long truncated;
        double seconds;

        @Override
        public String toString() {
            return String.format("packets=%d written=%d (IPv4 %d, IPv6 %d) skipped: non-IP %d, truncated %d; %.1f s",
                    packets, written, ipv4, ipv6, nonIp, truncated, seconds);
        }
    }
}
